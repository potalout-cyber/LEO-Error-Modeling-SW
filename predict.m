function [w1,w2,ca] = predict(altitude)

    % 원시 데이터 정의
    altitudex = [400, 600, 800, 1000, 1200, 1400];
    
    yw1 = [0.419, 0.487, 0.540, 0.582, 0.617, 0.647];
    yw2 = [0.642, 0.617, 0.595, 0.575, 0.556, 0.539];
    yca = [70.22, 66.07, 62.69, 59.82, 57.31, 55.09];
    
    % 보간/외삽 함수 생성 (선형 + 외삽 포함)
    f1 = @(xal) interp1(altitudex, yw1, xal, 'linear', 'extrap');
    f2 = @(xal) interp1(altitudex, yw2, xal, 'linear', 'extrap');
    f3 = @(xal) interp1(altitudex, yca, xal, 'linear', 'extrap');
    
    % 사용자가 원하는 x 값 입력
    alS = altitude;  % 여기에 원하는 x값들 입력
    
    % 함수 적용
    w1 = f1(alS);
    w2 = f2(alS);
    ca = f3(alS);

    assignin('base','wr',w1);
    assignin('base','wac',w2);
    assignin('base','ca',ca);
    
    % 결과 출력
    %disp('예측 결과:');
    %disp(table(alS', y1_pred', y2_pred', y3_pred', ...
    %    'VariableNames', {'orbit altitude(km)', 'Radial 가중치', 'Along,Cross 가중치', 'Cone Angle(deg)'}));
end