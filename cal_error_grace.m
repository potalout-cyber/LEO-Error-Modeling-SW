function [r, navPos, delta_O, O_error_RTN] = cal_error_grace(tle_filename, nasa_filename)
    
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


    %grace
    %leap_time = 27;
    nasa_unix_time = nasa_data{:, 1} + 946728000 - 18;

    nasaPos = nasa_data{:,7:9};
    nasaVel = nasa_data{:,13:15};

    
    [common_time, stk_idx, nasa_idx] = intersect(stk_unix_time, nasa_unix_time);
    
    % Error Vector
    delta_O = navPos(stk_idx,:) - nasaPos(nasa_idx,:);
    
    % 기준 위치와 속도
    r = nasaPos(nasa_idx,:);                % 위성 위치
    V = nasaVel(nasa_idx,:);                % 위성 속도 (ECEF)

    % 지구 자전 각속도 벡터 (rad/s)
    omega_e = [0, 0, 7.292115e-5];
    
    % omega_e x r 계산 (각 행마다)
    omega_cross_r = cross(repmat(omega_e, size(r,1), 1), r, 2);
    
    % 속도 벡터 v 계산
    v = V + omega_cross_r;
    
    % RTN 단위 벡터 계산
    e_R = normalize(r, 2);
    e_N = normalize(cross(r, v, 2), 2);
    e_T = normalize(cross(e_N, e_R, 2), 2);
    
    % 각 방향 성분 오차 계산
    O_error_R = dot(delta_O, e_R, 2);
    O_error_T = dot(delta_O, e_T, 2);
    O_error_N = dot(delta_O, e_N, 2);

    assignin('base','errorR',O_error_R);
    assignin('base','errorI',O_error_T);
    assignin('base','errorC',O_error_N);

    O_error_RTN = [O_error_R, O_error_T, O_error_N];
    
end