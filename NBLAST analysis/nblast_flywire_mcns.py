from __future__ import annotations

import argparse
import gc
import hashlib
import pickle
from pathlib import Path
import flybrains
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import navis
import pandas as pd

try:
	from scipy.optimize import linear_sum_assignment
except Exception:
	linear_sum_assignment = None

try:
	from scipy.cluster.hierarchy import fcluster, linkage, leaves_list
except Exception:
	fcluster = None
	linkage = None
	leaves_list = None


def _load_swc_dir(directory: Path, max_files: int | None = None) -> navis.NeuronList:
	if not directory.exists():
		raise FileNotFoundError(f"SWC directory not found: {directory}")

	swc_files = sorted(directory.glob("*.swc"))
	if max_files is not None:
		swc_files = swc_files[:max_files]
	if not swc_files:
		raise FileNotFoundError(f"No .swc files found in: {directory}")

	neurons = []
	for swc_file in swc_files:
		neuron = navis.read_swc(str(swc_file))
		try:
			neuron.id = int(swc_file.stem)
		except ValueError:
			neuron.id = swc_file.stem
		neurons.append(neuron)

	return navis.NeuronList(neurons)


def _sort_score_matrix(scores):
	row_order = scores.mean(axis=1).sort_values(ascending=False).index
	col_order = scores.mean(axis=0).sort_values(ascending=False).index
	return scores.loc[row_order, col_order]


def _safe_name(text: str) -> str:
	return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in str(text))


def _style_plot_axis(ax) -> None:
	ax.set_xlim(200000, 600000)
	ax.set_ylim(200000, 80000)
	ax.grid(False)
	for spine in ax.spines.values():
		spine.set_visible(False)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_axis_off()


def _find_one_to_one_pairs(scores: pd.DataFrame) -> pd.DataFrame:
	if scores.empty:
		return pd.DataFrame(columns=["flywire_id", "mcns_id", "score"])

	if linear_sum_assignment is not None:
		# Hungarian solves min-cost assignment; negate to maximize NBLAST score sum.
		cost = -scores.to_numpy()
		row_ind, col_ind = linear_sum_assignment(cost)
		pairs = pd.DataFrame(
			{
				"flywire_id": [scores.index[i] for i in row_ind],
				"mcns_id": [scores.columns[j] for j in col_ind],
				"score": [float(scores.iat[i, j]) for i, j in zip(row_ind, col_ind)],
			}
		)
		return pairs.sort_values("score", ascending=False).reset_index(drop=True)

	# Fallback: greedy one-to-one matching by highest remaining score.
	stacked = (
		scores.stack()
		.rename("score")
		.reset_index()
		.rename(columns={"level_0": "flywire_id", "level_1": "mcns_id"})
		.sort_values("score", ascending=False)
	)
	used_flywire = set()
	used_mcns = set()
	rows = []
	for row in stacked.itertuples(index=False):
		if row.flywire_id in used_flywire or row.mcns_id in used_mcns:
			continue
		rows.append({"flywire_id": row.flywire_id, "mcns_id": row.mcns_id, "score": float(row.score)})
		used_flywire.add(row.flywire_id)
		used_mcns.add(row.mcns_id)
		if len(used_flywire) >= min(scores.shape[0], scores.shape[1]):
			break
	return pd.DataFrame(rows)


def _validate_bijection(pairs: pd.DataFrame, scores: pd.DataFrame) -> None:
	n_flywire, n_mcns = scores.shape
	if n_flywire != n_mcns:
		raise ValueError(
			"Bijection requires equal counts, but got "
			f"FlyWire={n_flywire}, MCNS={n_mcns}."
		)

	if len(pairs) != n_flywire:
		raise ValueError(
			"Failed to build complete bijection: "
			f"expected {n_flywire} pairs, got {len(pairs)}."
		)

	if pairs["flywire_id"].nunique() != n_flywire:
		raise ValueError("Bijection check failed: duplicate FlyWire IDs in matched pairs.")
	if pairs["mcns_id"].nunique() != n_mcns:
		raise ValueError("Bijection check failed: duplicate MCNS IDs in matched pairs.")


def _cluster_and_group_scores(scores: pd.DataFrame, n_clusters: int):
	if scores.empty:
		empty = pd.DataFrame()
		return empty, empty, empty, empty

	if linkage is not None and leaves_list is not None and fcluster is not None and len(scores) > 1 and len(scores.columns) > 1:
		row_link = linkage(scores.to_numpy(), method="average", metric="euclidean")
		col_link = linkage(scores.to_numpy().T, method="average", metric="euclidean")
		row_order = [scores.index[i] for i in leaves_list(row_link)]
		col_order = [scores.columns[i] for i in leaves_list(col_link)]

		k_rows = max(1, min(n_clusters, len(scores.index)))
		k_cols = max(1, min(n_clusters, len(scores.columns)))
		row_labels_raw = fcluster(row_link, t=k_rows, criterion="maxclust")
		col_labels_raw = fcluster(col_link, t=k_cols, criterion="maxclust")
	else:
		# Fallback when scipy clustering is unavailable.
		row_order = scores.mean(axis=1).sort_values(ascending=False).index.tolist()
		col_order = scores.mean(axis=0).sort_values(ascending=False).index.tolist()
		row_labels_raw = [1] * len(scores.index)
		col_labels_raw = [1] * len(scores.columns)

	clustered_scores = scores.loc[row_order, col_order]

	flywire_labels = pd.DataFrame(
		{
			"flywire_id": scores.index,
			"cluster": [int(x) for x in row_labels_raw],
		}
	).sort_values(["cluster", "flywire_id"])
	mcns_labels = pd.DataFrame(
		{
			"mcns_id": scores.columns,
			"cluster": [int(x) for x in col_labels_raw],
		}
	).sort_values(["cluster", "mcns_id"])

	row_cluster_map = dict(zip(scores.index, flywire_labels.set_index("flywire_id")["cluster"]))
	col_cluster_map = dict(zip(scores.columns, mcns_labels.set_index("mcns_id")["cluster"]))

	grouped = scores.copy()
	grouped.index = [row_cluster_map[idx] for idx in grouped.index]
	grouped.columns = [col_cluster_map[idx] for idx in grouped.columns]
	grouped_scores = grouped.groupby(level=0).mean().groupby(level=0, axis=1).mean()
	grouped_scores.index = [f"flywire_cluster_{c}" for c in grouped_scores.index]
	grouped_scores.columns = [f"mcns_cluster_{c}" for c in grouped_scores.columns]

	return clustered_scores, flywire_labels, mcns_labels, grouped_scores


def _plot_correspondence_pairs(
	flywire_targets,
	mcns_targets,
	pairs: pd.DataFrame,
	output_dir: Path,
) -> None:
	plots_dir = output_dir / "pair_overlays"
	plots_dir.mkdir(parents=True, exist_ok=True)
	brain_mesh = flybrains.JRCFIB2022M.mesh_brain

	for idx, row in enumerate(pairs.itertuples(index=False), start=1):
		flywire_id = row.flywire_id
		mcns_id = row.mcns_id
		score = float(row.score)

		fig, ax = plt.subplots(1, 1, figsize=(8, 8))
		navis.plot2d(
			[flywire_targets.idx[flywire_id]*1000, mcns_targets.idx[mcns_id]*1000, brain_mesh],
			view=("x", "-z"),
			alpha=0.6,
			method="2d",
			ax=ax,
		)
		ax.set_title(f"Pair {idx}: FW={flywire_id} | MCNS={mcns_id} | score={score:.3f}")
		ax.set_aspect("equal")
		_style_plot_axis(ax)
		fig.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.03)
		plot_path = plots_dir / f"pair_{idx:03d}_{_safe_name(flywire_id)}__{_safe_name(mcns_id)}.png"
		fig.savefig(plot_path, dpi=600)
		fig.clf()
		plt.close(fig)
		del fig, ax
		if idx % 25 == 0:
			gc.collect()


def _plot_top5_panels_per_neuron(
	scores: pd.DataFrame,
	flywire_targets,
	mcns_targets,
	output_dir: Path,
) -> None:
	flywire_dir = output_dir / "top5_overlays_by_flywire"
	mcns_dir = output_dir / "top5_overlays_by_mcns"
	flywire_dir.mkdir(parents=True, exist_ok=True)
	mcns_dir.mkdir(parents=True, exist_ok=True)
	brain_mesh = flybrains.JRCFIB2022M.mesh_brain

	# For each FlyWire neuron: overlay with top-5 MCNS neurons.
	for idx, flywire_id in enumerate(scores.index, start=1):
		top5 = scores.loc[flywire_id].sort_values(ascending=False).head(5)
		fig_fw, axes_fw = plt.subplots(1, 5, figsize=(25, 5))
		for i, (mcns_id, score) in enumerate(top5.items()):
			ax = axes_fw[i]
			navis.plot2d(
				[
					flywire_targets.idx[flywire_id] * 1000,
					mcns_targets.idx[mcns_id] * 1000,
					brain_mesh,
				],
				view=("x", "-z"),
				alpha=0.6,
				method="2d",
				ax=ax,
			)
			ax.set_title(f"MCNS {mcns_id}\n{float(score):.3f}")
			ax.set_aspect("equal")
			_style_plot_axis(ax)

		for j in range(len(top5), 5):
			axes_fw[j].axis("off")

		fig_fw.suptitle(f"FlyWire {flywire_id} top-5 MCNS overlaps", y=0.98)
		fig_fw.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.04, wspace=0.08)
		fig_fw.savefig(flywire_dir / f"flywire_{_safe_name(flywire_id)}_top5.png", dpi=300)
		fig_fw.clf()
		plt.close(fig_fw)
		del fig_fw, axes_fw
		if idx % 20 == 0:
			plt.close("all")
			gc.collect()

	# For each MCNS neuron: overlay with top-5 FlyWire neurons.
	for idx, mcns_id in enumerate(scores.columns, start=1):
		top5 = scores.loc[:, mcns_id].sort_values(ascending=False).head(5)
		fig_mcns, axes_mcns = plt.subplots(1, 5, figsize=(25, 5))
		for i, (flywire_id, score) in enumerate(top5.items()):
			ax = axes_mcns[i]
			navis.plot2d(
				[
					mcns_targets.idx[mcns_id] * 1000,
					flywire_targets.idx[flywire_id] * 1000,
					brain_mesh,
				],
				view=("x", "-z"),
				alpha=0.6,
				method="2d",
				ax=ax,
			)
			ax.set_title(f"FlyWire {flywire_id}\n{float(score):.3f}")
			ax.set_aspect("equal")
			_style_plot_axis(ax)

		for j in range(len(top5), 5):
			axes_mcns[j].axis("off")

		fig_mcns.suptitle(f"MCNS {mcns_id} top-5 FlyWire overlaps", y=0.98)
		fig_mcns.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.04, wspace=0.08)
		fig_mcns.savefig(mcns_dir / f"mcns_{_safe_name(mcns_id)}_top5.png", dpi=300)
		fig_mcns.clf()
		plt.close(fig_mcns)
		del fig_mcns, axes_mcns
		if idx % 20 == 0:
			plt.close("all")
			gc.collect()

	plt.close("all")
	gc.collect()


def _dotprops_cache_path(
	cache_dir: Path,
	dataset_name: str,
	neurons: navis.NeuronList,
	k: int,
	resample: float,
	test: bool,
) -> Path:
	ids = ",".join(str(n.id) for n in neurons)
	fingerprint = hashlib.md5(ids.encode("utf-8")).hexdigest()[:10]
	resample_token = str(resample).replace(".", "p")
	test_token = "_test" if test else ""
	filename = f"{dataset_name}_k{k}_res{resample_token}_{fingerprint}{test_token}.pkl"
	return cache_dir / filename


def _load_or_make_dotprops(
	neurons: navis.NeuronList,
	dataset_name: str,
	cache_dir: Path,
	k: int,
	resample: float,
	test: bool,
):
	cache_path = _dotprops_cache_path(cache_dir, dataset_name, neurons, k, resample, test)
	if cache_path.exists():
		try:
			with cache_path.open("rb") as f:
				cached = pickle.load(f)
			print(f"Loaded cached dotprops for {dataset_name}: {cache_path}")
			return cached
		except Exception as exc:
			print(f"Failed to load dotprops cache {cache_path}: {exc}. Recomputing.")

	print(f"Building dotprops for {dataset_name} (k={k}, resample={resample})")
	dotprops = navis.make_dotprops(neurons, k=k, resample=resample)
	cache_path.parent.mkdir(parents=True, exist_ok=True)
	with cache_path.open("wb") as f:
		pickle.dump(dotprops, f, protocol=pickle.HIGHEST_PROTOCOL)
	print(f"Saved dotprops cache for {dataset_name}: {cache_path}")
	return dotprops


def run_nblast(
	swc_cache_dir: str = "./swc_cache",
	dotprops_cache_dir: str = "./swc_cache/dotprops",
	output_dir: str = "./outputs",
	dotprops_k: int = 5,
	dotprops_resample: float = 1000,
	nblast_scores: str = "mean",
	n_clusters: int = 6,
	test: bool = False,
) -> None:
	base = Path(swc_cache_dir)
	flywire_dir = base / "flywire"
	mcns_dir = base / "mcns"

	max_files = 3 if test else None
	flywire = _load_swc_dir(flywire_dir, max_files=max_files)
	mcns = _load_swc_dir(mcns_dir, max_files=max_files)

	if test:
		print("Test mode enabled: using first 3 SWC files from each dataset.")

	print(f"Loaded {len(flywire)} FlyWire SWCs from {flywire_dir}")
	print(f"Loaded {len(mcns)} MCNS SWCs from {mcns_dir}")

	dp_cache_path = Path(dotprops_cache_dir)
	flywire_dp = _load_or_make_dotprops(
		flywire,
		dataset_name="flywire",
		cache_dir=dp_cache_path,
		k=dotprops_k,
		resample=dotprops_resample/4,
		test=test,
	)
	mcns_dp = _load_or_make_dotprops(
		mcns,
		dataset_name="mcns",
		cache_dir=dp_cache_path,
		k=dotprops_k,
		resample=dotprops_resample/8,
		test=test,
	)
	flywire_xf = navis.xform_brain(flywire_dp, source="FLYWIRE", target="JRCFIB2022M")
	mcns_xf = navis.xform_brain(mcns_dp, source="JRCFIB2022Mraw", target="JRCFIB2022M")
	scores = navis.nblast(flywire_xf/1000, mcns_xf/1000, scores=nblast_scores, progress=True)
	sorted_scores = _sort_score_matrix(scores)
	pairs = _find_one_to_one_pairs(scores)
	_validate_bijection(pairs, scores)
	clustered_scores, flywire_labels, mcns_labels, grouped_scores = _cluster_and_group_scores(
		scores,
		n_clusters=n_clusters,
	)

	out = Path(output_dir)
	out.mkdir(parents=True, exist_ok=True)

	raw_path = out / "flywire_vs_mcns_scores.csv"
	sorted_path = out / "flywire_vs_mcns_sorted_scores.csv"
	pairs_path = out / "flywire_vs_mcns_one_to_one_pairs.csv"
	clustered_path = out / "flywire_vs_mcns_clustered_scores.csv"
	flywire_clusters_path = out / "flywire_cluster_labels.csv"
	mcns_clusters_path = out / "mcns_cluster_labels.csv"
	grouped_path = out / "flywire_vs_mcns_grouped_scores.csv"
	scores.to_csv(raw_path)
	sorted_scores.to_csv(sorted_path)
	pairs.to_csv(pairs_path, index=False)
	clustered_scores.to_csv(clustered_path)
	flywire_labels.to_csv(flywire_clusters_path, index=False)
	mcns_labels.to_csv(mcns_clusters_path, index=False)
	grouped_scores.to_csv(grouped_path)

	_plot_correspondence_pairs(flywire_xf / 1000, mcns_xf / 1000, pairs, out)
	_plot_top5_panels_per_neuron(scores, flywire_xf / 1000, mcns_xf / 1000, out)

	print(f"Saved raw score matrix: {raw_path}")
	print(f"Saved sorted score matrix: {sorted_path}")
	print(f"Saved one-to-one pairs: {pairs_path}")
	print(f"Saved clustered score matrix: {clustered_path}")
	print(f"Saved FlyWire cluster labels: {flywire_clusters_path}")
	print(f"Saved MCNS cluster labels: {mcns_clusters_path}")
	print(f"Saved grouped cluster score matrix: {grouped_path}")
	print(f"Saved pair overlay plots to: {out / 'pair_overlays'}")
	print(f"Saved FlyWire top-5 panel plots to: {out / 'top5_overlays_by_flywire'}")
	print(f"Saved MCNS top-5 panel plots to: {out / 'top5_overlays_by_mcns'}")


def _build_argparser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description=(
			"Load SWC files from ./swc_cache/flywire and ./swc_cache/mcns, "
			"run NBLAST, and save score matrices."
		)
	)
	parser.add_argument("--swc-cache-dir", default="./swc_cache", help="Base SWC cache directory")
	parser.add_argument(
		"--dotprops-cache-dir",
		default="./swc_cache/dotprops",
		help="Directory for cached dotprops files",
	)
	parser.add_argument("--output-dir", default="./outputs", help="Directory to save output CSV files")
	parser.add_argument("--dotprops-k", type=int, default=5, help="k-neighbors for dotprops")
	parser.add_argument(
		"--dotprops-resample",
		type=float,
		default=1000,
		help="Resample distance for dotprops generation",
	)
	parser.add_argument(
		"--nblast-scores",
		default="mean",
		choices=["forward", "mean", "min", "max", "both"],
		help="NBLAST scoring mode",
	)
	parser.add_argument(
		"--n-clusters",
		type=int,
		default=6,
		help="Number of clusters used to group FlyWire/MCNS neurons",
	)
	parser.add_argument(
		"--test",
		action="store_true",
		help="Run a quick test using only the first 3 SWC files from each dataset",
	)
	return parser


if __name__ == "__main__":
	args = _build_argparser().parse_args()
	run_nblast(
		swc_cache_dir=args.swc_cache_dir,
		dotprops_cache_dir=args.dotprops_cache_dir,
		output_dir=args.output_dir,
		dotprops_k=args.dotprops_k,
		dotprops_resample=args.dotprops_resample,
		nblast_scores=args.nblast_scores,
		n_clusters=args.n_clusters,
		test=args.test,
	)
