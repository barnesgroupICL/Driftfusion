function property = import_single_property(property_in, T, possible_headers, start_row, end_row)
error_checker = zeros(1, length(possible_headers));

for j = 1:length(possible_headers)
    try
        temp_param = T{: , possible_headers(j)}';
        property = temp_param(start_row:end_row);
    catch
        error_checker(j) = 1;
    end
end

if all(error_checker)
    %message = ['No column headings match:', possible_headers, ', using default in PC.'];
    %warning(char(message)')
    property = ones(1, end_row - start_row + 1)*property_in(1);
end
end