function [gregorian] = mjd2gregorian(mjd)

%     % Convert MJD to Julian Date (JD)
%     jd = mjd + 2400000.5;
%     
%     a = jd + 0.5;
% 
%     
%     % Step 4: Convert to Gregorian date using algorithm
%     if a < 2299161
%         c = a + 1524;
%     else
%         
%         b = floor((a - 1867216.25) / 36524.25);
%         c = a + b - floor(b/4) + 1525;
%     end
% 
%     d = floor((c-122.1)/365.25);
%     e = floor(365.25*d);
%     f = floor((c-e)/30.6001);
%     D = c - e - floor(30.6001*f) + (jd + 0.5 - a);
%     M = f - 1 - 12*floor(f/14);
%     Y = d - 4715 - floor((7+M)/10);
% 
%     % Display the result
%     fprintf('Gregorian Date: %d-%02d-%02.0f \n', Y, M, D);

    % Step 1: Convert MJD to Julian Date (JD)
    jd = mjd + 2400000.5;
    
    % Step 2: Add 0.5 to JD for conversion (accounts for the fact that JD starts at noon)
    jd = jd + 0.5;

    % Step 3: Break JD into integer and fractional parts
    Z = floor(jd); % Integer part
    F = jd - Z;    % Fractional part
    
    % Step 4: Convert to Gregorian date using algorithm
    if Z < 2299161
        A = Z;
    else
        alpha = floor((Z - 1867216.25) ./ 36524.25);
        A = Z + 1 + alpha - floor(alpha ./ 4);
    end

    B = A + 1524;
    C = floor((B - 122.1) ./ 365.25);
    D = floor(365.25 .* C);
    E = floor((B - D) ./ 30.6001);
    
    % Step 5: Day calculation
    day = B - D - floor(30.6001 .* E) + F;

    % Step 6: Month and year calculation
    if E < 14
        month = E - 1;
    else
        month = E - 13;
    end

    if month > 2
        year = C - 4716;
    else
        year = C - 4715;
    end

    % Step 7: Extract hour, minute, second from the fractional day
    fractional_day = F .* 24;
    hour = floor(fractional_day);
    minute = floor((fractional_day - hour) .* 60);
    second = (fractional_day - hour - minute ./ 60) .* 3600;
    
    gregorian = [year, month, floor(day), hour, minute, second];

    % Display the result
%     fprintf('Gregorian Date: %d-%02d-%02.0f %02d:%02d:%06.3f\n', year, month, floor(day), hour, minute, second);


end

