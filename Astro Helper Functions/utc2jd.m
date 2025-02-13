function [JD] = utc2jd(year,month,day,hour,minute,second)

    % Outputs
    % sidereal_angle: [degrees]
    
    % If the month is January or February, adjust year and month
    if month <= 2
        year = year - 1;
        month = month + 12;
    end
    
    % Calculate A and B
    A = floor(year / 100);
    B = 2 - A + floor(A / 4);
    
    % Calculate the Julian Day
    JD = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;
    
    % Compute fractional day
    frac_minutes = minute + floor(second/60) + mod(second,60)/60;
    frac_hours = hour + floor(frac_minutes/60) + mod(frac_minutes,60)/60;
    fractional_day = floor(frac_hours/24) + mod(frac_hours,24)/24;
    JD = JD + fractional_day;

end