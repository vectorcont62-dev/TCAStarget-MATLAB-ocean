%% ========================================
%% getRadiusIndex.m
%% cell_info_initial의 radius를 Rad.r 그리드 인덱스로 변환
%% ========================================

function [idx, radius_quantized, error] = getRadiusIndex(radius_pos, Rad, varargin)
%
% 함수 설명:
% cell_info_initial에서의 radius 값을 Rad.r의 이산 그리드에서의 인덱스로 변환
%
% Input:
%   radius_pos [double]: cell_info_initial{i,j}(1)의 값 (m)
%   Rad [struct]: setRadius 함수의 출력
%   (선택) show_info [logical]: 디버그 정보 출력 여부 (기본값: false)
%
% Output:
%   idx [int]: Rad.r 배열의 인덱스 (1 ≤ idx ≤ length(Rad.r))
%   radius_quantized [double]: 실제 양자화된 위치값 = Rad.r(idx)
%   error [double]: 근사 위치와 양자화된 위치의 차이

    % 입력 파싱
    p = inputParser;
    addParameter(p, 'show_info', false, @islogical);
    parse(p, varargin{:});
    show_info = p.Results.show_info;
    
    % ========== 직접 계산 방식 (빠름) ==========
    % 수식: idx = round((radius - r_min) / dr) + 1
    
    r_min = Rad.r(1);
    dr = Rad.dr;
    N = length(Rad.r);
    
    % 인덱스 계산
    idx_raw = round((radius_pos - r_min) / dr) + 1;
    
    % 경계 체크 (중요!)
    if idx_raw < 1
        idx = 1;
        warning_flag = true;
    elseif idx_raw > N
        idx = N;
        warning_flag = true;
    else
        idx = idx_raw;
        warning_flag = false;
    end
    
    % 실제 양자화된 위치값
    radius_quantized = Rad.r(idx);
    
    % 오차
    error = abs(radius_quantized - radius_pos);
    
    % 디버그 정보 (선택사항)
    if show_info
        layer_name = determineLayer(radius_pos, Rad);
        fprintf('[getRadiusIndex] %s: r_input=%.3e m → idx=%d, r_grid=%.3e m, error=%.3e m\n', ...
            layer_name, radius_pos, idx, radius_quantized, error);
        if warning_flag
            fprintf('  ⚠️ 경계 조정됨 (raw_idx=%d)\n', idx_raw);
        end
    end
end

%% 보조 함수: 위치로부터 층 판단
function layer_name = determineLayer(radius_pos, Rad)
    if radius_pos < Rad.r_o1
        layer_name = 'O1';
    elseif radius_pos < Rad.r_n1
        layer_name = 'N1';
    elseif radius_pos < Rad.r_o2
        layer_name = 'O2';
    elseif radius_pos < Rad.r_ctn
        layer_name = 'CTN';
    else
        layer_name = 'BOX/OTHER';
    end
end
