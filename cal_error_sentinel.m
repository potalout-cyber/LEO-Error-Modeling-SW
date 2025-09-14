function [r, navPos, delta_O] = cal_error_sentinel(tle_filename, nasa_filename)
    
    % 파일 옵션 설정: 모든 열을 문자열로 먼저 불러오기
    opts = detectImportOptions(tle_filename);
    opts = setvartype(opts, 'char');
    tle_data = readtable(tle_filename, opts);
    
    % UTC 기준 시간 문자열을 datetime으로 변환
    datetime_format = 'dd MMM yyyy HH:mm:ss.SSS'; 
    time_strings = tle_data{:,1};
    t_utc = datetime(time_strings, 'InputFormat', datetime_format, ...
                     'Locale', 'en_US', 'TimeZone', 'UTC');
    
    % Unix Time 변환
    stk_unix_time = posixtime(t_utc);
    
    % Position 데이터 (x, y, z) 추출 후 숫자 변환
    navPos = str2double(tle_data{:, 2:4});
    navVel = str2double(tle_data{:, 5:7});

    
    % NASA 데이터
    nasa_data = readtable(nasa_filename);

    %sentinel
    nasa_unix_time = nasa_data{:, 1};
    nasaPos = nasa_data{:,2:4};
    
    [common_time, stk_idx, nasa_idx] = intersect(stk_unix_time, nasa_unix_time);
    
    % Error Vector
    delta_O = navPos(stk_idx,:) - nasaPos(nasa_idx,:);
    
    % 기준 위치와 속도
    r = nasaPos(nasa_idx,:);                % 위성 위치
    
end