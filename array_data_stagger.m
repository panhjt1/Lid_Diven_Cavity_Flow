classdef array_data_stagger
    %ARRAY_DATA Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ncell_X
        ncell_Y

        X
        Y

        U_stg
        V_stg

        T
        P
        
        Sc % 源项
        Sp % 线性化源项，可不管，但是稳定性会更差
        
        rho
        lam_vis
        lam_cond
    end
    properties(Hidden) % 修正量隐藏不然调用的时候很乱
        Pprime
    end
    properties(Hidden,Constant)
        % 一些方便使用的函数句柄直接做成常量放到这
        mat_e = @(mat) [mat(:, 2:end), NaN * mat(:, 1)];
        mat_w = @(mat) [NaN * mat(:, 1), mat(:, 1:end-1)];
        mat_s = @(mat) [mat(2:end, :); NaN * mat(1, :)];
        mat_n = @(mat) [NaN * mat(1, :); mat(1:end-1, :)];

        scheme_CD     = @(Peclet_d) 1 - 0.5*Peclet_d;
        scheme_FUD    = @(Peclet_d) Peclet_d.^0; % 这个写法是为了Pe矩阵里有NaN时矩阵里边对应体现出NaN，后边的零矩阵同理
        scheme_Hybrid = @(Peclet_d) max(Peclet_d.^0 - 1, 1 - 0.5*Peclet_d);
        scheme_Exp    = @(Peclet_d) Peclet_d / (exp(Peclet_d) - 1);
        scheme_Power  = @(Peclet_d) max(Peclet_d.^0 - 1, (1 - 0.1*Peclet_d).^5);
        % TODO：无法这样写，得写成函数
        interp_n = @(rho) (rho .* obj.mat_n(dys) + obj.mat_n(rho) .* dyn) ./ (obj.mat_n(dys) + dyn);% 从原位网格的格点向北侧边界插值
        interp_s = @(rho) (rho .* obj.mat_s(dyn) + obj.mat_s(rho) .* dys) ./ (obj.mat_s(dyn) + dys);
        interp_e = @(rho) (rho .* obj.mat_e(dyw) + obj.mat_e(rho) .* dye) ./ (obj.mat_e(dyw) + dye);
        interp_w = @(rho) (rho .* obj.mat_w(dye) + obj.mat_w(rho) .* dyw) ./ (obj.mat_w(dye) + dyw);
    end
    properties(Hidden, Dependent)
        dxe % 原位网格E方向至界面距离，以下相似
        dxw
        dyn
        dys
        XU % 交错U网格格点位置
        YV % 交错V网格格点位置
        nnode_X
        nnode_Y
        Vol
        VolU
        VolV
    end

    methods(Static)
        function [sparse_coef, b] = coef_mat_to_sparse(aP, aE, aW, aN, aS, b)
            [m, n] = size(aP);
            % 展开所有矩阵
            aP = aP.';
            aP = aP(:);

            aE = aE.';
            aE = aE(:);

            aW = aW.';
            aW = aW(:);

            aN = aN.';
            aN = aN(:);

            aS = aS.';
            aS = aS(:);
            
            b = b.';
            b = b(:);
            
            sparse_coef = spdiags([-aS, -aE, aP, -aW, -aN], [-n -1:1 n], m*n, m*n).';
        end

    end

    methods
        function obj = array_data_stagger(ncell_X, ncell_Y, X_domain, Y_domain)
            %ARRAY_DATA Construct an instance of this class
            %   Detailed explanation goes here
            obj.ncell_X = ncell_X;
            obj.ncell_Y = ncell_Y;

            obj.lam_vis  = ones(obj.nnode_X, obj.nnode_Y);
            obj.lam_cond = ones(obj.nnode_X, obj.nnode_Y);

            dx = X_domain / obj.ncell_X;
            dy = Y_domain / obj.ncell_Y;
            X = [0,dx/2:dx:X_domain-dx/2, X_domain];
            Y = [Y_domain,Y_domain-dy/2:-dy:dy/2, 0];
            [obj.X, obj.Y] = meshgrid(X, Y);
            
            obj.P = zeros(obj.nnode_X, obj.nnode_Y);
            obj.U_stg = zeros(obj.nnode_X, obj.nnode_Y);
            obj.V_stg = zeros(obj.nnode_X, obj.nnode_Y);
            obj.rho = ones(obj.nnode_X, obj.nnode_Y);
            obj.Sc   = zeros(obj.nnode_X, obj.nnode_Y);
            obj.Sp   = zeros(obj.nnode_X, obj.nnode_Y);
        end
        %%
        function nnode_X = get.nnode_X(obj)
            nnode_X = obj.ncell_X + 2;
        end
        function nnode_Y = get.nnode_Y(obj)
            nnode_Y = obj.ncell_Y + 2;
        end
        function dxe = get.dxe(obj)
            dx = abs(diff(obj.X,1, 2)); % 原位网格控制容积的dx，同样用指向的原节点命名
            dx = [NaN*dx(:,1), dx];
            dxe = [0 * dx(:, end), dx(:, 3:end-1)/2, dx(:, end), 0 * dx(:, end)];% 原位网格dxe
        end
        function dxw = get.dxw(obj)
            dx = abs(diff(obj.X,1, 2)); % 原位网格控制容积的dx，同样用指向的原节点命名
            dx = [NaN*dx(:,1), dx];
            dxw = [0 * dx(:, end), dx(:, 2), dx(:, 3:end-1)/2, 0 * dx(:, end)];
        end
        function dyn = get.dyn(obj)
            dy = abs(diff(obj.Y, 1, 1));
            dy = [dy; NaN*dy(1,:)];
            dyn = [0 * dy(1, :); dy(1, :); dy(2:end-2, :)/2; 0 * dy(end-1, :)];
        end
        function dys = get.dys(obj) 
            dy = abs(diff(obj.Y, 1, 1));
            dy = [dy; NaN*dy(1,:)];
            dys = [0 * dy(1, :); dy(2:end-2, :)/2; dy(end-1, :); 0 * dy(end-1, :)];
        end
        function XU = get.XU(obj)
            XU = obj.X - obj.dxw; % 交错U网格格点位置
            XU(:,1) = NaN;
        end
        function YV = get.YV(obj)
            YV = obj.Y - obj.dys; % 交错V网格格点位置
            YV(end, :) = NaN;
        end
        function Vol = get.Vol(obj)
            Vol = (obj.dxw + obj.dxe) .* (obj.dyn + obj.dys);
        end
        function VolU = get.VolU(obj)
            VolU_dx = obj.dxw + obj.mat_w(obj.dxe);
            VolU_dx(:, 3) = VolU_dx(:, 3) + obj.dxw(:, 2);
            VolU_dx(:, end-1) = VolU_dx(:, end-1) + obj.dxe(:, end-1);
            VolU = VolU_dx .* (obj.dyn + obj.dys);
            VolU(:, [2, end]) = 0;
            VolU([1, end], :) = 0;
            VolU(:, 1) = NaN;
        end
        function VolV = get.VolV(obj)
            VolV_dy = obj.dys + obj.mat_s(obj.dyn);
            VolV_dy(2, :) = VolV_dy(2, :) + obj.dyn(2,:);
            VolV_dy(end-2, :) = VolV_dy(end-2, :) + obj.dys(end-1, :);
            VolV = (obj.dxw + obj.dxe) .* VolV_dy;
            VolV([1,end-1],:) = 0;
            VolV(:, [1,end]) = 0;
            VolV(end, :) = NaN;
        end
        %%
        function phi_n = f_n(obj, phi)
            phi_n = (phi .* obj.mat_n(obj.dys) + obj.mat_n(phi) .* obj.dyn) ./ (obj.mat_n(obj.dys) + obj.dyn);% 从原位网格的格点向北侧边界插值
        end
        
        function phi_s = f_s(obj, phi)
            phi_s = (phi .* obj.mat_s(obj.dyn) + obj.mat_s(phi) .* obj.dys) ./ (obj.mat_s(obj.dyn) + obj.dys);
        end
        
        function phi_e = f_e(obj, phi)
            phi_e = (phi .* obj.mat_e(obj.dxw) + obj.mat_e(phi) .* obj.dxe) ./ (obj.mat_e(obj.dxw) + obj.dxe);
            
        end
        function phi_w = f_w(obj, phi)
            phi_w = (phi .* obj.mat_w(obj.dxe) + obj.mat_w(phi) .* obj.dxw) ./ (obj.mat_w(obj.dxe) + obj.dxw);
        end





        function [D_U_X, D_U_Y, F_U_X, F_U_Y] = flux_Ugrid(obj, GAMMA)
            % 计算扩散通量和对流通量D、F
            dxe = obj.dxe;
            dxw = obj.dxw;
            dyn = obj.dyn;
            dys = obj.dys;

            GAMMA_X = GAMMA; % 每个交错U容积的右边界值
            GAMMA_X(:, 2) = GAMMA(:, 1); % 每个交错U容积的右边界值
            GAMMA_X(:, end-1) = GAMMA_X(:, end);
            GAMMA_X(:, 1) = NaN * GAMMA(:, 1); % 每个交错U容积的右边界值
            D_U_X = GAMMA_X .* (dyn + dys) ./ (dxw + dxe);% 检查ok
            D_U_X(:,end) = NaN;
            F_U_X = obj.rho .* (obj.U_stg .* dxe + obj.mat_e(obj.U_stg) .* dxw) ./ (dxe + dxw)...
                .* (dyn + dys); % 简化，部分密度直接取节点值不插值，速度线性插值
            % FUX两侧处理
            F_U_X(:, 2) = obj.rho(:, 1) .* obj.U_stg(:, 2) .* (dyn(:,2) + dys(:,2));
            F_U_X(:, end-1) = obj.rho(:, end) .* obj.U_stg(:, end) .* (dyn(:,end-1) + dys(:,end-1));% 检查ok

            dxe(:, 2) = dxe(:, 2) + dxw(:, 2); % 左侧处理
            dxw(:,end-1) = dxe(:, end-1) + dxw(:, end-1); % 右侧处理
            D_U_Y = obj.mat_w(dxe) ./ (dyn ./ obj.mat_w(GAMMA) + obj.mat_n(dys) ./ obj.mat_w(obj.mat_n(GAMMA))) ...
                + dxw            ./ (dyn ./ GAMMA            + obj.mat_n(dys) ./ obj.mat_n(GAMMA)); % 检查ok
%             interp_n = @(rho) (rho .* obj.mat_n(dys) + obj.mat_n(rho) .* dyn) ./ (obj.mat_n(dys) + dyn); % 与N节点之间线性插值
            F_U_Y = obj.f_n(obj.rho) .* obj.mat_n(obj.V_stg) .* dxw ...
                + obj.f_n(obj.mat_w(obj.rho)) .* obj.mat_w(obj.mat_n(obj.V_stg)) .* obj.mat_w(dxe); % 检查ok
            % FUY左右两侧处理,似乎不需要？

        end
        
        function [D_V_X, D_V_Y, F_V_X, F_V_Y] = flux_Vgrid(obj, GAMMA)
            % V grid
            dxe = obj.dxe;
            dxw = obj.dxw;
            dyn = obj.dyn;
            dys = obj.dys;

            GAMMA_Y = GAMMA; % V容积上边界值
            GAMMA_Y(2, :) = GAMMA_Y(1,:);
            GAMMA_Y(end-1, :) = GAMMA_Y(end,:);
            GAMMA_Y(end, :) = NaN * GAMMA_Y(end,:);
            D_V_Y = GAMMA_Y .* (dxe + dxw) ./ (dyn + dys);
            D_V_Y(1,:) = NaN;
            F_V_Y = obj.rho .* (obj.V_stg .* obj.mat_n(dys) + obj.mat_n(obj.V_stg) .* dyn) ./ (dyn + obj.mat_n(dys))...
                .* (dxe + dxw);
            % FVY上下处理
            F_V_Y(2, :) = obj.rho(1, :) .* obj.V_stg(1, :) .* (dxe(2, :) + dxw(2, :));
            F_V_Y(end-1,:)= obj.rho(end, :) .* obj.V_stg(end-1,:)  .* (dxe(end-1, :) + dxw(end-1, :));
            F_V_Y(end, :) = NaN;

            dys(2, :) = dyn(2, :) + dys(2, :);
            dyn(end-1, :) = dyn(end-1, :) + dys(end-1, :);
            D_V_X = dys          ./ (dxe ./ GAMMA            + obj.mat_e(dxw) ./ obj.mat_e(GAMMA))...
                + obj.mat_s(dyn) ./ (dxe ./ obj.mat_s(GAMMA) + obj.mat_e(dxw) ./ obj.mat_s(obj.mat_e(GAMMA)));
%             interp_e = @(rho) (rho .* obj.mat_e(dxw) + obj.mat_e(rho) .* dxe) ./ (obj.mat_e(dxw) + dxe);
            F_V_X = obj.f_e(obj.rho)            .* obj.mat_e(obj.U_stg) .* dys...
                + obj.f_e(obj.mat_s(obj.rho)) .* obj.mat_s(obj.mat_e(obj.U_stg)) .* obj.mat_s(dyn);
        end

        function [aP, aE, aW, aN, aS, b] = coef(obj, D_X, D_Y, F_X, F_Y, dt, phi_name, scheme)
            dxe = obj.dxe;
            dxw = obj.dxw;
            dyn = obj.dyn;
            dys = obj.dys;

            Pe_X = F_X ./ D_X;
            Pe_Y = F_Y ./ D_Y;
            aE = scheme(abs(Pe_X)) .* D_X + max(-F_X, 0);
            aE(isnan(aE)) = 0;
            % aW根据aE计算后调整得来
            aW = F_X + aE;
            aW(:, 2:end) = aW(:, 1:end-1);
            aW(:, 1) = 0;
            aW(isnan(aW)) = 0;
            Pe_Y = F_Y ./ D_Y;
            aN = scheme(abs(Pe_Y)) .* D_Y + max(-F_Y, 0);
            aN(isnan(aN)) = 0;
            % aS根据aN计算后调整得来
            aS = F_Y + aN;
            aS(1:end-1, :) = aS(2:end, :);
            aS(end, :) = 0;
            aS(isnan(aS)) = 0;
            if strcmp(phi_name, "U")
                phi = obj.U_stg;
                Vol = obj.VolU;
                DP = obj.P - obj.mat_w(obj.P);
                DP(:, 2)    = 0;
                DP(:, end)  = 0;
                DP(:, 3)    = DP(:, 3)     .* (dxe(:, 2)     + dxw(:, 2)     + dxw(:, 3)) ./ (dxe(:, 2)    + dxw(:, 3));%边缘单元格外插
                DP(:, end-1)= DP(:, end-1) .* (dxe(:, end-1) + dxw(:, end-1) + dxe(:, end-2)) ./ (dxw(:, end-1) + dxe(:, end-2));
                DX          = obj.X - obj.mat_w(obj.X); 
                DX(:, 2)    = 0;
                DX(:, end)  = 0;
                DX(:, 3)    = dxe(:, 2)     + dxw(:, 2)     + dxw(:, 3);
                DX(:, end-1)= dxe(:, end-1) + dxw(:, end-1) + dxe(:, end-2);
                DP_DXi = DP ./ DX;
            elseif strcmp(phi_name, "V")
                phi = obj.V_stg;
                Vol = obj.VolV;
                DP = obj.P - obj.mat_s(obj.P);
                DP(1, :)    = 0;
                DP(end-1, :)= 0;
                DP(2, :)    = DP(2, :)     .* (dyn(2, :) + dys(2, :) + dyn(3, :)) ./ (dys(2, :) + dyn(3, :));
                DP(end-2, :)= DP(end-2, :) .* (dyn(end-1, :) + dys(end-1, :) + dys(end-2, :)) ./ (dyn(end-1, :) + dys(end-2, :));
                DY          = obj.Y - obj.mat_s(obj.Y); 
                DY(1, :)    = 0;
                DY(end-1, :)= 0;
                DY(2, :)    = dyn(2, :) + dys(2, :) + dyn(3, :);
                DY(end-2, :)= dyn(end-1, :) + dys(end-1, :) + dys(end-2, :);
                DP_DXi = DP ./ DY;
            elseif strcmp(phi_name, "T")
                phi = obj.T;
                Vol = obj.Vol;
                DP_DXi = 0;
                % TODO：未完成
            else % 错误
                error("未识别到变量");
            end

            aP0 = obj.rho .* Vol / dt;
            aP = aE + aW + aS + aN + aP0 - obj.Sp .* Vol;                  
            aP(isnan(aP)) = 0;

            b = (obj.Sc - DP_DXi) .* Vol + aP0 .* phi; % p_term为压力相关项，对于温度等变量为0，只有动量方程中非零

            % 边界值处理
            [aP, aE, aW, aN, aS, b] = obj.bound(aP, aE, aW, aN, aS, b, phi_name);
        end

        function [aP_P, aE_P, aW_P, aN_P, aS_P, b_P] = coef_Pprime(obj, aP_U, aP_V)
            % 输入动量方程的系数，输出求p'的修正方程系数
            dxe = obj.dxe;
            dxw = obj.dxw;
            dyn = obj.dyn;
            dys = obj.dys;

            de = (dyn + dys) ./ obj.mat_e(aP_U);
            dw = (dyn + dys) ./ aP_U;
            dn = (dxe + dxw) ./ obj.mat_n(aP_V);
            ds = (dxe + dxw) ./ aP_V;
            % 界面上的rho进行插值
            rho_n = obj.f_n(obj.rho);
            rho_s = obj.f_s(obj.rho);
            rho_e = obj.f_e(obj.rho);
            rho_w = obj.f_w(obj.rho);
            aE_P = de .* (dyn + dys) .* rho_e;% 这个地方看起来可以跳过求de，暂时先放着方便自己读
            aW_P = dw .* (dyn + dys) .* rho_w;
            aN_P = dn .* (dxe + dxw) .* rho_n;
            aS_P = ds .* (dxe + dxw) .* rho_s;
            b_P = (rho_w .* obj.U_stg - rho_e .* obj.mat_e(obj.U_stg)) .* (dyn + dys) + ...
                  (rho_s .* obj.V_stg - rho_n .* obj.mat_n(obj.V_stg)) .* (dxe + dxw); 
            % 以上求b假设常物性，所以认为rhoP-rhoP0=0
            % 此处书上rhoP-rhoP0为连续性方程导出，似乎对于压力基都应该为0？
            % 而且b收敛时最终归零，所以应该只影响数值稳定性不影响结果
            
            % 处理边界条件
            [~, aE_P, aW_P, aN_P, aS_P, b_P] = obj.bound([], aE_P, aW_P, aN_P, aS_P, b_P, "Pprime"); % 这里限制之后再aE+aW+aS+aN才是正确的aP，因为实现第一类边界条件的时候四边的系数被人工修改过，直接相加会有问题
            aP_P = aE_P + aW_P + aN_P + aS_P;

        end

        function [aP, aE, aW, aN, aS, b] = bound(obj, aP, aE, aW, aN, aS, b, phi_name)
            % 把所有边界都按第一类边界处理，改一下系数。如果有其他类型边界用附加源项法
            tmp = cat(3, aE, aW, aN, aS);
            tmp(:, [1, end], :) = 0;
            tmp([1, end], :, :) = 0;
            if strcmp(phi_name, "U")
                tmp(:, 2, :) = 0;
                tmp(:, 1, :) = NaN;
                aP(:, [2, end]) = 1;
                aP([1, end], :) = 1;
                aP(:, 1, :) = NaN;
                b(:, [2, end]) = obj.U_stg(:, [2, end]);
                b([1, end], :) = obj.U_stg([1, end], :);
                b(:, 1, :)=NaN;
            elseif strcmp(phi_name, "V")
                tmp(end-1, :, :) = 0;
                tmp(end, :, :) = NaN;
                aP(:, [1, end]) = 1;
                aP([1, end-1],:) = 1;
                aP(end, :, :) = NaN;
                b(:, [1, end]) = obj.V_stg(:, [1, end]);
                b([1, end-1],:) = obj.V_stg([1, end-1], :);
            elseif strcmp(phi_name, "Pprime")
                % 动量修正，四边不参与运算，
                tmp(:, end-1, 1) = 0;% aE
                tmp(:,     2, 2) = 0;% aW
                tmp(2,     :, 3) = 0;% aN
                tmp(end-1, :, 4) = 0;% aS
            elseif strcmp(phi_name, "T")
                aP(:, [1, end]) = 1;
                aP([1, end],:) = 1;
                b(:, [1, end]) = obj.T(:, [1, end]);
                b([1, end],:) = obj.T([1, end],:);
            end
            aE = tmp(:,:,1);
            aW = tmp(:,:,2);
            aN = tmp(:,:,3);
            aS = tmp(:,:,4);
        end
        function result_mat = solve_eqn(obj, aP, aE, aW, aN, aS, b, phi_name)
            % 把东南西北系数转换成矩阵方程并求解
            % 对于U、V网格，把方程里边不需要的行列去掉再求解，结果里边重新加回来保证作图方便
            % obj.coef_mat_to_sparse
            tmp = cat(3, aP, aE, aW, aN, aS, b);
            if strcmp(phi_name, "U")
                tmp(:, 1, :) = [];
            elseif strcmp(phi_name, "V")
                tmp(end, :, :) = [];
            elseif strcmp(phi_name, "Pprime")
                tmp([1, end], :, :) = [];
                tmp(:, [1, end], :) = [];
            else % 错误
                error("未识别到变量");
            end

            aP = tmp(:,:,1);
            aE = tmp(:,:,2);
            aW = tmp(:,:,3);
            aN = tmp(:,:,4);
            aS = tmp(:,:,5);
            b  = tmp(:,:,6);
            [h,w] = size(aP);
            [sparse_coef, b] = obj.coef_mat_to_sparse(aP, aE, aW, aN, aS, b);
            result_mat = reshape(sparse_coef\b, [w h]).';
            
            if strcmp(phi_name, "U")
                result_mat = [result_mat(:,1)*NaN, result_mat];
            elseif strcmp(phi_name, "V")
                result_mat = [result_mat; result_mat(end,:)*NaN];
            elseif strcmp(phi_name, "Pprime")
%                 result_mat = [0 * result_mat(:,1), result_mat, 0 * result_mat(:,end)];
%                 result_mat = [0 * result_mat(1,:); result_mat; 0 * result_mat(end,:)];
                result_mat = [nan * result_mat(:,1), result_mat, nan * result_mat(:,end)];
                result_mat = [nan * result_mat(1,:); result_mat; nan * result_mat(end,:)];
            else % 错误
                error("未识别到变量");
            end

        end

        function [Uprime, Vprime] = vel_correct(obj, aP_U, aP_V)
            dxe = obj.dxe;
            dxw = obj.dxw;
            dyn = obj.dyn;
            dys = obj.dys;

            de = (dyn + dys) ./ aP_U;
            ds = (dxe + dxw) ./ aP_V;
            
            de(:,  [2,end])=0;% 没有完全想明白这里是不是应该这样置0，似乎是
            ds([1,end-1],:)=0;

            Uprime = (obj.mat_w(obj.Pprime) - obj.Pprime) .* de; % 此处的p'之差在边缘处应不应该外插？需要！
            Uprime(:, 3)    = Uprime(:, 3)    .* (dxe(:, 2) + dxw(:, 2) + dxw(:, 3)) ./ (dxe(:, 2) + dxw(:, 3));%边缘单元格外插
            Uprime(:, end-1)= Uprime(:, end-1).* (dxe(:, end-1) + dxw(:, end-1) + dxe(:, end-2)) ./ (dxw(:, end-1) + dxe(:, end-2));
            % 速度第一类边界条件，四边修正值为0
            Uprime(:,  2)   = 0;
            Uprime(:,  end) = 0;
            Uprime(1,  :)   = 0;
            Uprime(end,:)   = 0;
            Uprime(:,1)     = nan;
            Vprime = (obj.mat_s(obj.Pprime) - obj.Pprime) .* ds;
            Vprime(2, :) =Vprime(2, :) .* (dyn(2, :) + dys(2, :) + dyn(3, :)) ./ (dys(2, :) + dyn(3, :));
            Vprime(end-2, :)= Vprime(end-2, :) .* (dyn(end-1, :) + dys(end-1, :) + dys(end-2, :)) ./ (dyn(end-1, :) + dys(end-2, :));
            % 速度第一类边界条件，四边修正值为0
            Vprime(:,  1)   = 0;
            Vprime(:,  end) = 0;
            Vprime(1,  :)   = 0;
            Vprime(end-1,:) = 0;
            Vprime(end,:)   = nan;
        end


        function [UU, VV] = vel_interp(obj)
            % 把交错网格的速度插值回原位网格上
            UU = griddata(...
                obj.XU(~isnan(obj.XU)), obj.Y(~isnan(obj.XU)), obj.U_stg(~isnan(obj.XU)),...
                obj.X, obj.Y,'linear');
            VV = griddata(...
                obj.X(~isnan(obj.YV)), obj.YV(~isnan(obj.YV)), obj.V_stg(~isnan(obj.YV)),...
                obj.X, obj.Y,'linear');
        end
    end
end
