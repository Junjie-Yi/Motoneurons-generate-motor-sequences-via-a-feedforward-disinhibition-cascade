from __future__ import annotations

import argparse
import gc
import hashlib
import pickle
from pathlib import Path
import flybrains
import matplotlib.pyplot as plt
import navis


def _safe_name(text: str) -> str:
	return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in text)


def _style_plot_axis(ax, view) -> None:
	ax.set_xlim(200000, 600000)
	if view == ("x", "-z"):
		ax.set_ylim(200000, 80000)
	elif view == ("x", "-y"):
		ax.set_ylim(400000, 200000)
	ax.grid(False)
	for spine in ax.spines.values():
		spine.set_visible(False)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_axis_off()


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


def _load_lm_query(nrrd_file: Path, query_k: int, threshold: int, thin: bool):
	def _get_n_points(dp) -> int:
		n_points = getattr(dp, "n_points", None)
		if n_points is not None:
			return int(n_points)
		points = getattr(dp, "points", None)
		if points is not None:
			return int(len(points))
		return 0

	original_threshold = threshold
	current_threshold = threshold
	max_retries = 10

	query = navis.read_nrrd(
		str(nrrd_file),
		output="dotprops",
		threshold=current_threshold,
		thin=thin,
		k=query_k,
	)
	n_points = _get_n_points(query)

	retries = 0
	while n_points > 50000 and retries < max_retries:
		current_threshold = current_threshold + original_threshold
		retries += 1
		print(
			f"Query {nrrd_file.stem}: n_points={n_points} > 50000, "
			f"reloading with threshold={current_threshold}"
		)
		query = navis.read_nrrd(
			str(nrrd_file),
			output="dotprops",
			threshold=current_threshold,
			thin=thin,
			k=query_k,
		)
		n_points = _get_n_points(query)

	if n_points > 50000:
		print(
			f"Warning: query {nrrd_file.stem} still has n_points={n_points} "
			f"after {max_retries} retries; using threshold={current_threshold}."
		)

	query.id = nrrd_file.stem
	query.name = nrrd_file.stem
	return query


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
	dotprops = navis.make_dotprops(neurons, k=k, parallel=True, resample=resample)
	cache_path.parent.mkdir(parents=True, exist_ok=True)
	with cache_path.open("wb") as f:
		pickle.dump(dotprops, f, protocol=pickle.HIGHEST_PROTOCOL)
	print(f"Saved dotprops cache for {dataset_name}: {cache_path}")
	return dotprops


def _plot_top10(query, target_dotprops, top10_scores, out_dir: Path) -> None:
	cm = 1 / 2.54
	fig = plt.figure(figsize=(15 * cm, 6 * cm))
	
	import matplotlib.gridspec as gridspec
	outer_gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.3) # Approximate relative spacing for 0.5cm
	
	axes_flat = []
	for outer_idx in range(2):
		inner_gs = gridspec.GridSpecFromSubplotSpec(2, 5, subplot_spec=outer_gs[outer_idx], hspace=0.1, wspace=0.017, height_ratios=[0.6, 1.0])
		for row_idx in range(2):
			for col_idx in range(5):
				axes_flat.append(fig.add_subplot(inner_gs[row_idx, col_idx]))

	views = [("x", "-z"), ("x", "-y")]

	# i is 0-9
	for i, (target_id, score) in enumerate(top10_scores.items()):
		outer_row = i // 5
		col = i % 5
		
		for v_idx, view in enumerate(views):
			# Index in axes_flat: outer_row * 10 + v_idx * 5 + col
			ax_idx = outer_row * 10 + v_idx * 5 + col
			ax = axes_flat[ax_idx]
			
			target = target_dotprops.idx[target_id]
			navis.plot2d(
				[query * 1000, target, flybrains.JRCFIB2022M.mesh_brain],
				view=view,
				alpha=0.6,
				linewidth=0.2,
				method="2d",
				ax=ax,
			)
			
			if v_idx == 0: # Only plot title on the ("x", "-z") top view
				ax.set_title(f"{target_id} | {float(score):.3f}", fontsize=5, pad=2)
			
			ax.set_aspect("equal")
			_style_plot_axis(ax, view)

	# Turn off unused axes
	for j in range(len(top10_scores), 10):
		outer_row = j // 5
		col = j % 5
		for v_idx in range(2):
			axes_flat[outer_row * 10 + v_idx * 5 + col].axis("off")

	fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.05)
	fig.savefig(out_dir / "top10_overlay.png", dpi=600)
	plt.close(fig)


def _run_one_dataset(
	query,
	query_name: str,
	dataset_name: str,
	target_dotprops,
	output_root: Path,
	nblast_scores: str,
) -> None:
	scores = navis.nblast(query, target_dotprops/1000, scores=nblast_scores, progress=True)
	ranked = scores.loc[query_name].sort_values(ascending=False)
	top10 = ranked.head(10)

	out_dir = output_root / _safe_name(query_name) / dataset_name
	out_dir.mkdir(parents=True, exist_ok=True)

	scores.to_csv(out_dir / "nblast_scores_full.csv")
	top10.to_csv(out_dir / "top10_scores.csv", header=["score"])
	_plot_top10(query, target_dotprops, top10, out_dir)

	print(f"Saved {dataset_name} results for {query_name} -> {out_dir}")


def run_batch(
	mn_dir: str = "./mn",
	swc_cache_dir: str = "./swc_cache",
	dotprops_cache_dir: str = "./swc_cache/dotprops",
	output_dir: str = "./outputs/lm_vs_em",
	query_k: int = 20,
	query_threshold: int = 200,
	query_thin: bool = True,
	target_dotprops_k: int = 5,
	target_resample_um: float = 1000/8,
	nblast_scores: str = "forward",
	test: bool = False,
	mn12v_scale: float | None = None,
) -> None:
	mn_path = Path(mn_dir)
	lm_files = sorted(mn_path.glob("*.nrrd"))
	if test:
		lm_files = lm_files[:3]
	if not lm_files:
		raise FileNotFoundError(f"No .nrrd files found in: {mn_path}")

	max_swc_files = 3 if test else None
	flywire_swc = _load_swc_dir(Path(swc_cache_dir) / "flywire", max_files=max_swc_files)
	manc_swc = _load_swc_dir(Path(swc_cache_dir) / "manc", max_files=max_swc_files)

	if test:
		print("Test mode enabled: using first 3 LM NRRD files and first 3 SWC files per EM dataset.")

	print(f"Loaded {len(flywire_swc)} FlyWire SWCs")
	print(f"Loaded {len(manc_swc)} MCNS SWCs")

	dp_cache_path = Path(dotprops_cache_dir)
	flywire_dp = _load_or_make_dotprops(
		flywire_swc,
		dataset_name="flywire",
		cache_dir=dp_cache_path,
		k=target_dotprops_k,
		resample=target_resample_um,
		test=test,
	)
	del flywire_swc
	gc.collect()

	manc_dp = _load_or_make_dotprops(
		manc_swc,
		dataset_name="manc",
		cache_dir=dp_cache_path,
		k=target_dotprops_k,
		resample=target_resample_um,
		test=test,
	)
	del manc_swc
	gc.collect()

	flywire_xf = navis.xform_brain(flywire_dp, source="FLYWIRE", target="JRCFIB2022M")
	del flywire_dp
	gc.collect()
	
	manc_xf = navis.xform_brain(manc_dp, source="JRCFIB2022Mraw", target="JRCFIB2022M")
	del manc_dp
	gc.collect()

	output_root = Path(output_dir)
	output_root.mkdir(parents=True, exist_ok=True)

	for nrrd_file in lm_files:
		query = _load_lm_query(nrrd_file, query_k=query_k, threshold=query_threshold, thin=query_thin)
		query_name = str(query.id)
		if query_name == "MN12V" and mn12v_scale is not None:
			print(f"Rescaling MN12V by {mn12v_scale}")
			query = query * mn12v_scale
		query_xf = navis.xform_brain(query, source="JRC2018U", target="JRCFIB2022Mum")
		del query
		gc.collect()

		print(f"Running LM query: {query_name}")
		_run_one_dataset(
			query=query_xf,
			query_name=query_name,
			dataset_name="flywire",
			target_dotprops=flywire_xf,
			output_root=output_root,
			nblast_scores=nblast_scores,
		)
		_run_one_dataset(
			query=query_xf,
			query_name=query_name,
			dataset_name="mcns",
			target_dotprops=manc_xf,
			output_root=output_root,
			nblast_scores=nblast_scores,
		)
		del query_xf
		gc.collect()


def _build_argparser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description=(
			"Load LM NRRD files from ./mn, compare each LM to FlyWire and MCNS SWC datasets "
			"from ./swc_cache, and save top-10 scores and figures per LM."
		)
	)
	parser.add_argument("--mn-dir", default="./mn", help="Directory containing LM .nrrd files")
	parser.add_argument("--swc-cache-dir", default="./swc_cache", help="SWC cache base directory")
	parser.add_argument(
		"--dotprops-cache-dir",
		default="./swc_cache/dotprops",
		help="Directory for cached dotprops files",
	)
	parser.add_argument("--output-dir", default="./outputs/lm_vs_em", help="Output root directory")
	parser.add_argument("--query-k", type=int, default=20, help="k-neighbors for LM dotprops")
	parser.add_argument("--query-threshold", type=int, default=300, help="NRRD threshold for LM query")
	parser.add_argument("--query-thin", action="store_true", default=True, help="Enable thinning")
	parser.add_argument("--no-query-thin", action="store_false", dest="query_thin", help="Disable thinning")
	parser.add_argument("--target-dotprops-k", type=int, default=5, help="k-neighbors for EM dotprops")
	parser.add_argument("--target-resample-um", type=float, default=1000/4, help="Resample size for EM dotprops")
	parser.add_argument(
		"--nblast-scores",
		default="mean",
		choices=["forward", "mean", "min", "max", "both"],
		help="NBLAST scoring mode",
	)
	parser.add_argument(
		"--test",
		action="store_true",
		help="Run a quick test with first 3 LM files and first 3 SWC files per EM dataset",
	)
	parser.add_argument(
		"--mn12v-scale",
		type=float,
		default=None,
		help="Factor to rescale MN12V LM data if provided (e.g. 1.0 or other values)",
	)
	return parser


if __name__ == "__main__":
	args = _build_argparser().parse_args()
	run_batch(
		mn_dir=args.mn_dir,
		swc_cache_dir=args.swc_cache_dir,
		dotprops_cache_dir=args.dotprops_cache_dir,
		output_dir=args.output_dir,
		query_k=args.query_k,
		query_threshold=args.query_threshold,
		query_thin=args.query_thin,
		target_dotprops_k=args.target_dotprops_k,
		target_resample_um=args.target_resample_um,
		nblast_scores=args.nblast_scores,
		test=args.test,
		mn12v_scale=args.mn12v_scale,
	)
