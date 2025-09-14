function fetch_and_save_TLE(NORADID, startDate, endDate)
    % 사용자 로그인 정보
    username = 'ID@naver.com';
    password = 'PW';
    cookiesFile = 'cookies.txt';

    % 날짜 문자열 변환
    startStr = datestr(startDate, 'yyyy-mm-dd');
    endStr   = datestr(endDate, 'yyyy-mm-dd');

    % 로그인
    loginCmd = sprintf('curl -X POST -d "identity=%s&password=%s" --cookie-jar %s https://www.space-track.org/ajaxauth/login', ...
                        username, password, cookiesFile);
    [status, loginOut] = system(loginCmd);
    if status ~= 0
        error("❌ 로그인 실패: %s", loginOut);
    end

    % TLE 다운로드
    tleURL = sprintf(['https://www.space-track.org/basicspacedata/query/class/tle/NORAD_CAT_ID/%s/EPOCH/%s--%s/orderby/EPOCH%%20asc/format/tle'], ...
                      NORADID, startStr, endStr);
    tleCmd = sprintf('curl --cookie %s "%s" -o tle_data.txt', cookiesFile, tleURL);
    [status, tleOut] = system(tleCmd);
    if status ~= 0
        error("❌ TLE 다운로드 실패: %s", tleOut);
    end

    % TLE 유효성 검사
    lines = readlines('tle_data.txt');
    lines = lines(~cellfun(@isempty, cellstr(lines)));
    if mod(length(lines), 2) ~= 0
        error('❌ TLE 형식 오류: 짝수 줄 아님');
    end
    fprintf("✅ TLE 다운로드 성공 (%d줄)\n", length(lines));

    % 상태벡터 계산
    computeStateVectors(NORADID, startDate, endDate, 'tle_data.txt');
end

function computeStateVectors(NORADID, startTime, endTime, tleFile)
    sampleInterval = 1;
    outputCSV = 'tle_data/satellite_state_vectors.csv';

    try
        sc = satelliteScenario(startTime, endTime, sampleInterval);
        sat = satellite(sc, tleFile, 'Name', NORADID, 'OrbitPropagator', 'sgp4');
    catch
        error('❌ satelliteScenario 생성 실패: TLE 파일 오류 또는 파싱 실패');
    end

    timeVec = startTime:seconds(sampleInterval):endTime;
    nTimes = numel(timeVec);
    pos_ecef = zeros(nTimes, 3);
    vel_ecef = zeros(nTimes, 3);

    for i = 1:nTimes
        [r, v] = states(sat, timeVec(i));
        jd = juliandate(timeVec(i));
        [r_ecef, v_ecef] = teme2ecef(r, v, jd);
        pos_ecef(i, :) = r_ecef;
        vel_ecef(i, :) = v_ecef;
    end

    timeStr = cellstr(datestr(timeVec', 'dd mmm yyyy HH:MM:SS.FFF'));
    resultTable = table(timeStr, ...
        pos_ecef(:,1), pos_ecef(:,2), pos_ecef(:,3), ...
        vel_ecef(:,1), vel_ecef(:,2), vel_ecef(:,3), ...
        'VariableNames', {'Time_UTCG', 'x (m)', 'y (m)', 'z (m)', ...
                          'vx (m/sec)', 'vy (m/sec)', 'vz (m/sec)'});

    writetable(resultTable, outputCSV);
    fprintf("✔️ 상태벡터 저장 완료: %s\n", outputCSV);
end

function [r_ecef, v_ecef] = teme2ecef(r_teme, v_teme, jd)
    gmst = siderealTime(jd);
    R = [cos(gmst), sin(gmst), 0;
        -sin(gmst), cos(gmst), 0;
         0, 0, 1];
    r_ecef = (R * r_teme(:))';
    v_ecef = (R * v_teme(:))';
end

function theta = siderealTime(jd)
    T = (jd - 2451545.0) / 36525.0;
    theta_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T + ...
                0.093104*T^2 - 6.2e-6*T^3;
    theta_sec = mod(theta_sec, 86400);
    theta = theta_sec * (pi / 43200);  % rad
end
