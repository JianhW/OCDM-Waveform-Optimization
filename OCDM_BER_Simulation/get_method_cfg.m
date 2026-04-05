function cfg = get_method_cfg(method)
%GET_METHOD_CFG 统一管理不同算法方法的接口配置
% method 可以是方法名，也可以是已构造好的配置结构体

if isstruct(method)
    cfg = method;
    return;
end

if isstring(method)
    method = char(method);
end

switch method
    case 'OCDM_raw'
        cfg.key = 'OCDM_raw';
        cfg.displayName = 'OCDM';
        cfg.txTransform = 'IDFnT';
        cfg.rxUseY = true;
        cfg.rxPostIFFT = true;

    case 'OFDM_raw'
        cfg.key = 'OFDM_raw';
        cfg.displayName = 'OFDM';
        cfg.txTransform = 'IFFT';
        cfg.rxUseY = false;
        cfg.rxPostIFFT = false;

    otherwise
        error('Unknown method: %s', method);
end
end