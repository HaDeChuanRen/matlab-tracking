classdef CACFARdet
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        n_train
        n_guard
        CFAR_win
        num_cells
        dim_win
    end

    methods
        function obj = CACFARdet(in_tra, in_gua)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.n_train = in_tra;
            obj.dim_win = length(obj.n_train);
            obj.n_guard = [in_gua; zeros(obj.dim_win - length(in_gua), 1)];

%             for d_idx = 1 : obj.dim_win
%                 if obj.n_guard(d_idx) > obj.n_train(d_idx)
%                     obj.n_guard(d_idx) = obj.n_train(d_idx);
%                 end
%             end
            if obj.dim_win == 1
                obj.CFAR_win = ones(2 * (obj.n_train + obj.n_guard) + 1, 1);
                obj.CFAR_win(obj.n_train + (1 : (2 * obj.n_guard + 1))) = 0;
                obj.num_cells = (2 * (obj.n_train + obj.n_guard) + 1) -...
                    (2 * obj.n_guard + 1);
            elseif obj.dim_win == 2
                obj.CFAR_win = ones(2 * (obj.n_train(1) + obj.n_guard(1))...
                    + 1, 2 * (obj.n_train(2) + obj.n_guard(2)) + 1);
                obj.CFAR_win(obj.n_train(1) + (1 : (2 * obj.n_guard(1) + 1)),...
                    obj.n_train(2) + (1 : (2 * obj.n_guard(2) + 1))) = 0;
                obj.num_cells = (2 * (obj.n_train(1) + obj.n_guard(1)) + 1)...
                    * (2 * (obj.n_train(2) + obj.n_guard(2)) + 1) -...
                    (2 * obj.n_guard(1) + 1) * (2 * obj.n_guard(2) + 1);
            else
                error("3D以上的没做哈哈哈哈~");
            end
        end

        function mea_det = detection2D(obj, inMat, alpha_dB)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [Nx, My] = size(inMat);
            if obj.dim_win == 1
                obj.n_train = [obj.n_train; 0];
                obj.n_guard = [obj.n_guard; 0];
            end
            mat_pad = padarray(inMat, ...
                [obj.n_train(1) + obj.n_guard(1), ...
                obj.n_train(2) + obj.n_guard(2)], ...
                "symmetric", "both");
            conv_pad = conv2(mat_pad, obj.CFAR_win);
            sigma_mat = conv_pad(obj.n_train(1) + obj.n_guard(1) + (1 : Nx),...
                obj.n_train(2) + obj.n_guard(2) + (1 : My)) / obj.num_cells;
            snr_mat = 10 * log10(inMat ./ sigma_mat);
            mea_idx = find(snr_mat > alpha_dB);
            [mea_x, mea_y] = ind2sub([Nx, My], mea_idx);
            mea_det = [mea_x, mea_y];
        end
    end
end