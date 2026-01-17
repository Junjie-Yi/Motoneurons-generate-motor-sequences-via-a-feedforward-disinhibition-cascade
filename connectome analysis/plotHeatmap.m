warning off
[m,m_log2,data] = readAdj("matrix_Premotor_to_MN_min5.csv");
m_filtered = m;
m_filtered(sum(abs(m),2)<40,:) = [];
cg = clustergram(m_filtered,'Standardize',3,'DisplayRange',200,'Symmetric',true,'Colormap',redbluecmap,'Cluster',1,'linkage','median');
for ii = 1:9
    for jj  = 1:9
        ss(ii,jj) = cosineSimilarity(m_filtered(:,ii),m_filtered(:,jj));
    end
end
imagesc(ss);

function [m,m_norm,data] = readAdj(fileName)
    data = readtable(fileName,'ReadVariableNames',true,'ReadRowNames',true);
    m = data{:,:};
    m(abs(m)<5) = 0;
    m_sign = 2 * (m>-1) - 1;
    m_norm = log2(abs(m)+1).*m_sign;
end
function s = cosineSimilarity(a,b)
    s = sum((a.*b)) / (sqrt(sum(a.^2)) * sqrt(sum(b.^2)));
end