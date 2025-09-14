%gridsearch_URE.m

function [gridURE_maxs, gridmean_URE, gridmax_URE, gridstd_URE, pd] = get_grid_URE(r, navPos, alttitude, delta_O)

% 지표면 격자 생성
lat_range = -90:1:90;
lon_range = -180:1:180;
[lon_grid, lat_grid] = meshgrid(lon_range, lat_range);

% 격자 WGS84 타원체 기준 ECEF 좌표로 변환 (고도는 0으로 설정) , error나서 'meters'를 ()로 변경
[x_ecef, y_ecef, z_ecef] = geodetic2ecef(wgs84Ellipsoid('meters'), lat_grid, lon_grid, zeros(size(lat_grid)));


% 파라미터 설정
idx_size = size(r,1);           % 데이터 길이
[~, ~, y3_pred] = predict (alttitude); % 고도 기반 콘앵글 추정
cone_angle = y3_pred % y3_pred : 콘앵글
gridURE_maxs = zeros(idx_size, 1); % 결과 저장할 배열

% 계산 실행 
for i = 1:idx_size % idx만큼 실행(데이터 전체)
    
    sat_true = r(i,:)';     % i에서 위성 실제 위치 (nasa 데이터)
    sat_est = navPos(i,:)'; % i에서 위성 추정 위치 (tle->sgp4 데이터)
    error_vec = delta_O(i,:)'; % 시간 i에서 오차벡터 행렬

    % 커버리지 마스킹 함수 호출
    coverage_mask = coverage_masking(sat_true, x_ecef, y_ecef, z_ecef, cone_angle);
    
    % 격자점->위성 벡터
    griddx = sat_true(1) - x_ecef;
    griddy = sat_true(2) - y_ecef;
    griddz = sat_true(3) - z_ecef;
    
    % 벡터 크기 계산
    norm_vec = sqrt(griddx.^2 + griddy.^2 + griddz.^2);
    
    % 격자점에서 LoS 단위벡터 (LoS 벡터/크기)
    ux = griddx ./ norm_vec;
    uy = griddy ./ norm_vec;
    uz = griddz ./ norm_vec;
    
    % 각 격자점에서 LoS로 투영 URE = error_vec · LoS 단위벡터
    dot_product = error_vec(1) * ux + error_vec(2) * uy + error_vec(3) * uz;
    URE_grid = abs(dot_product);        % 절댓값
    URE_grid(~coverage_mask) = NaN;     % 커버리지가 아닌 경우 NaN

    % 평균 URE 계산 (NaN값 제외)
    gridURE_maxs(i) = max(URE_grid(~isnan(URE_grid)));
end

% 통계 및 시각화
gridmean_URE = mean(gridURE_maxs, 'omitnan');
gridmax_URE = max(gridURE_maxs, [], 'omitnan');
gridstd_URE = std(gridURE_maxs, 'omitnan');

fprintf("커버리지 max URE: %.3f m\n", gridmean_URE);
fprintf("최대: %.3f m\n", gridmax_URE);
fprintf("표준편차(σ) : %.3f m\n", gridstd_URE);

%figure;
%histogram(URE_means);
%title("커버리지 내 평균 URE 히스토그램");
%xlabel("URE [m]");
%ylabel("빈도");

% 반환값 선언
%assignin('base', 'mean_URE', mean_URE);
%assignin('base', 'std_URE', std_URE);    
%assignin('base', 'max_URE', max_URE);  
%assignin('base', 'URE_means', URE_means); 
%assignin('base', 'coverage_mask', coverage_mask);
assignin('base', 'gridURE_means', gridURE_maxs);

% 가우시안

% NaN 제거한 URE 데이터 추출
data = gridURE_maxs;

% 정규분포 피팅
pd = fitdist(data, 'Normal');  % 또는 'normal'

% x축 범위 설정 및 정규분포 PDF 계산
x = linspace(min(data), max(data), 200);
y = pdf(pd, x);  % 정규분포 pdf 계산

% 히스토그램 + 정규분포 곡선 시각화
figure;
histogram(data, 'Normalization', 'pdf'); hold on;
plot(x, y, 'LineWidth', 2);
title("Grid Calculation User Range Error");
xlabel("URE [m]");
ylabel("확률 밀도 (PDF)");
legend('히스토그램', '정규분포 피팅');

end