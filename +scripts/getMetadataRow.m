function [value, comment] = getMetadataRow(keywords, key)
idx = ismember(keywords(:, 1), key);
value = 0;
comment = 0;
if any(idx)
	value = keywords{idx, 2};
	comment = keywords{idx, 3};
end

end