function coverage_mask = coverage_masking(sat_pos, x_ecef, y_ecef, z_ecef, min_elvangle)
% 위성이 관측 가능한 격자 커버리지 계산

% 위성 방향 벡터 (각 격자점 → 위성)
vec_to_sat_x = sat_pos(1) - x_ecef;
vec_to_sat_y = sat_pos(2) - y_ecef;
vec_to_sat_z = sat_pos(3) - z_ecef;

% 지표점 → 위성 방향 벡터 계산 후 고도각(Elevation angle) 산출
dot_product = vec_to_sat_x .* x_ecef + vec_to_sat_y .* y_ecef + vec_to_sat_z .* z_ecef; % 대수 내적
norm_grid = sqrt(x_ecef.^2 + y_ecef.^2 + z_ecef.^2);                                    % ecef 벡터 크기
norm_vec = sqrt(vec_to_sat_x.^2 + vec_to_sat_y.^2 + vec_to_sat_z.^2);                   % sat 벡터 크기
cos_theta = dot_product ./ (norm_grid .* norm_vec);                                     % 기하 내적 변환으로 cos쎄타값 구하기
cos_theta = min(max(cos_theta, -1), 1); % 부동소수점 오차 방지 max, min 1, -1로 지정
elevation = asin(cos_theta) * (180/pi);     
% elevation = 90 - theta(지구 중심-위성 벡터와 사용자 시선 방향 벡터 사이의 각)
% cos(theta) = sin(90 - theta) = sin(elevation)
% arcsin(sin(elevation)) = elevation

% 최소 고도각 넘는 커버리지 추출->마스킹
coverage_mask = elevation > min_elvangle; % elevation 최소고도각보다 큰 경우 true

end