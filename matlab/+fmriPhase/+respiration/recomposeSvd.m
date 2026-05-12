function [volume, raw] = recomposeSvd(U, S, V, componentIdx, spatialSize)
%RECOMPOSESVD Recompose a volume from one SVD component.

componentMask = zeros(size(S));
componentMask(componentIdx, componentIdx) = 1;
raw = (S .* componentMask) * V';
volume = reshape(U * raw, [spatialSize, size(U, 2)]);

end
