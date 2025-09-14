function [gwoure_wl_all, max_oure, gwrms_oure]=get_worst_URE(r, navPos, delta_O)

XL_all = r;
Xerr_all = delta_O;
gwidx_size = size(r,1);           % 데이터 길이

% === OURE 계산 ===
gwoure_wl_all = zeros(gwidx_size, 1);

for i = 1:gwidx_size
    gwoure_val = calc_OURE_wl_precise(XL_all(i,:)', Xerr_all(i,:)');
    gwoure_wl_all(i) = gwoure_val;
end

assignin('base', 'URE_means', gwoure_wl_all);

figure;
plot(1:gwidx_size, gwoure_wl_all, 'b.-'); grid on;
xlabel('idx'); ylabel('OURE_{wl} (m)');
title('Worst-Case OURE over Time');

[max_oure, max_idx] = max(abs(gwoure_wl_all));
gwrms_oure = sqrt(mean(gwoure_wl_all.^2));

fprintf('\n\n🔹 Worst-case OURE_wl = %.3f m at index %d (UTC Time = %s)\n', max_oure, max_idx);
fprintf('🔸 RMS OURE_wl = %.3f m\n', gwrms_oure);

% NaN 제거한 URE 데이터 추출
gwoure_wl_all = gwoure_wl_all(~isnan(gwoure_wl_all));

% 정규분포 피팅
pd = fitdist(gwoure_wl_all, 'Normal');  % 또는 'normal'

% x축 범위 설정 및 정규분포 PDF 계산
x = linspace(min(gwoure_wl_all), max(gwoure_wl_all), 200);
y = pdf(pd, x);  % 정규분포 pdf 계산

% 히스토그램 + 정규분포 곡선 시각화
figure;
histogram(gwoure_wl_all, 'Normalization', 'pdf'); hold on;
plot(x, y, 'LineWidth', 2);
title("Worst Case User Range Error");
xlabel("URE [m]");
ylabel("확률 밀도 (PDF)");
legend('히스토그램', '정규분포 피팅');

% === URA 계산 (p = 1e-5) ===
%p = 1e-5;
%C_inv = @(p) sqrt(2) * erfinv(2 * p - 1);  % 정규분포 역함수 대체
%z_p = C_inv(1 - p/2);
%oura = rms_oure / z_p;

%fprintf('🔸 URA (%.1e integrity risk) = %.3f m\n', p, oura);

% === OURE 계산 함수 ===
function oure_wl = calc_OURE_wl_precise(XL, Xerr)
    f = 1 / 298.257222101;
    a = 6378137.0;
    e2 = 2*f - f^2;

    xL = XL(1); yL = XL(2); zL = XL(3);
    xerr = Xerr(1); yerr = Xerr(2); zerr = Xerr(3);

    beta2 = (xerr^2 + yerr^2)/(a^2 * zerr^2) + 1/(a^2 * (1 - f)^2);
    beta1 = (2 / (a^2 * zerr)) * (xerr * (xL - (xerr * zL / zerr)) + yerr * (yL - (yerr * zL / zerr)));
    beta0 = ((zerr * xL - xerr * zL)^2 + (zerr * yL - yerr * zL)^2) / (a^2 * zerr^2) - 1;
    D = beta1^2 - 4 * beta2 * beta0;

    if D >= 0
        dotprod = dot(XL, Xerr);
        oure_wl = sign(dotprod) * norm(Xerr);
    else
        max_cos_theta = -1;
        for phi_deg = -90:0.1:90
            phi = deg2rad(phi_deg);
            cos_phi = cos(phi);
            if abs(cos_phi) < 1e-6
                continue;
            end
            sin_phi = sin(phi);
            N_phi = a / sqrt(1 - e2 * sin_phi^2);

            % === 논문 식 (24) 기반 설명 추가 ===
            % Eq. (24): sin(lambda + arctan(b2/b1)) = b3 / sqrt(b1^2 + b2^2)
            % 여기서 lambda = asin(b3 / sqrt(b1^2 + b2^2)) - atan2(b2, b1)
            % 이는 위 식을 변형한 형태로, 논문과 동일한 위상 정렬 구조임
            b1 = yL;
            b2 = xL;
            b3 = a * sqrt(1 - e2 * sin_phi^2) / cos_phi - tan(phi) * zL;

            lambda = asin(b3 / sqrt(b1^2 + b2^2)) - atan2(b2, b1);

            cos_lambda = cos(lambda);
            sin_lambda = sin(lambda);
            xT = N_phi * cos_phi * cos_lambda;
            yT = N_phi * cos_phi * sin_lambda;
            zT = N_phi * (1 - e2) * sin_phi;
            XT = [xT; yT; zT];

            vec1 = (XL - XT) / norm(XL - XT);
            vec2 = Xerr / norm(Xerr);
            cos_theta = dot(vec1, vec2);
            if abs(cos_theta) > abs(max_cos_theta)
                max_cos_theta = cos_theta;
            end
        end
        oure_wl = norm(Xerr) * abs(max_cos_theta);
    end
end

end