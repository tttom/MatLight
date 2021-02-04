%
% [x_optim, err] = splitBregmanMinimize(x0, l1_function, l2_function, minimize_l2, l1_substitute_convergence_weight, nb_iterations, stop_criterion, want_statistics, ground_truth, progress_callback)
%
% Finds the x that minimizes ||l1_function(x)||_1 for ||l2_function(x)||_2^2<epsilon
% where epsilon converges to zero with increasing iteration.
%
% This implementation is based on the algorithm description in:
% Tom Goldstein and Stanley Osher "The Split Bregman Method for L1-Regularized Problems"
% SIAM J. Imaging Sci., 2(2), 323–343. https://doi.org/10.1137/080725891
%
%
% Inputs:
%     x0: a starting point for the solution. This can be a matrix or anything else that can be
%     passed to the l1_function and l2_function as an argument. If
%     ground_truth is specified it should be the same object and
%     subtractable to enable comparison.
%
%     l1_function: a function that takes a matrix or struct of the same size as
%         x0 as input and returns a 'shrinkable' object such as a column vector
%         of real values. Complex values should be split in two columns:
%         [real(x), imag(x)] to l1-minimize abs(x).
%         A 'shrinkable' object is anything to which the function or method
%         shrink(obj, l1Weight) can be applied, resulting in reduced shrinkable
%         object of the same dimensions but with values, u, that minimize
%         l1Weight*||joint(u)||_1+||u-obj||_2^2, where joint(u) is either the
%         identity or a function that joins parts of the input using the l2-norm
%         (isotropic shrinking, e.g. of a gradient vector length). When a
%         vector is specified, joint(u) := u, and when a matrix is specified,
%         joint(u) := sqrt(abs(sum(u, 2)).^2), i.e. rows are considered vectors.
%         Instead of a matrix, an instance of Shrinkable can be used as well.
%         This allows you to l1-minimize mixtures of scalar and vectors.
%
%     l2_function: returns a vector to be minimized in the l2-norm.
%
%     minimize_l2: A function with parameters (l1_target,l2_target_offset, l1SubstConvergenceWeight[, x0, nb_iterations(3)])
%             finds the u that minimizes l1SubstConvergenceWeight^2*||l1_function(u)-l1_target||_2^2 + ||l2_function(u)-l2_target_offset||_2^2
%             When an empty matrix [] is provided non-linear optimization
%             is performed automatically. The number of iterations can
%             be configured by the argument value nb_iterations(3).
%             An optional 4th input argument can be taken as an approximate
%             starting value for the optimization algorithm.
%             An optional 5th input argument specifies the desired number
%             of iterations as given by the user in nb_iterations(3), default 1.
%
%     l1_substitute_convergence_weight:
%             A scalar indicating the weight of the l1 substitute to make sure that l1-norm converges before the l1_substitute
%             Default: 0.10
%
%     nb_iterations: A vector with the number of outer and inner iterations to perform.
%                   A third value can indicate the iterations performed in
%                   the default minimize_l2.
%
%     stop_criterion(outerIt, l1_error, l2_error):
%                   a function that returns true when the outer iteration
%                   is to be stopped before nb_iterations(1) are reached. It
%                   can accept the input arguments outerIt, l1_error, and
%                   l2_error, where outerIt is the numer of outer iterations
%                   currently performed, and l1_error is the mean absolute
%                   error of l1_function and l2_error is root-mean-square error
%                   of the l2_function.
%                   When equal to the empty list, the iteration is stopped
%                   when the outer iteration reaches nb_iterations(1).
%                   The default is [].
%
%     ground_truth: (Optional) The ideal solution, when provided this is used to
%     track the evolution of the solution.
%
%     progress_callback: a call back function executed every iteration, taking the arguments:
%           iterationNb: a vector specifying the current outer to inner iteration numbers
%           totalNbIterations: a vector specifying the total number of iterations that will be performed per loop
%           restored: the current optimization result
%           l1_substitute: the substitute variable for the l1 function
%           l1_Bregman_constrainers: the Bregman variable to minimize the
%                   difference between the l1 substitute and the l1 function
%           l2_Bregman_constrainer: the Bregman variable to minimize the l2 function
%
% Contact Tom Vettenburg <tom.vettenburg@gmail.com> for more information.
% Don't leave this in if reusing this code or parts of it. Thanks.
%
% Outputs:
%     x_optim: The result of size identical to x0
%     err: a structure with several matrixes containing the l1 and l2
%     error, as well as its sum.
%
% Examples:
%     % find the peaks in 1D by minimizing the number of values that is non-zero
%     variableSize = [50, 1];  % A 1D problem
%     l1Weight = 1/0.05;  % how important is sparsity compared to best fit?
%     l1_subst_weight = 0.5;  % controls convergence rate
%     x0 = zeros(variableSize);
%     ground_truth = zeros(variableSize);
%     ground_truth([15, 25, 45]) = [1, 0.5, 0.75];
%     observation = ground_truth + 0.05*randn(variableSize);
%     l1_function = @(x) l1Weight*(x); % sparsity of the signal => spikes
%     l2_function = @(x) x-observation; % denoising, no deblurring
%     minimize_l2 = [];
%     nb_iterations = [50 10];  % The number of outer and inner iterations to do
%     [x_optim, err] = splitBregmanMinimize(x0, l1_function, l2_function, minimize_l2, l1_subst_weight, nb_iterations, [], true, ground_truth);
%     itNumbers = 1/nb_iterations(2)+[0:1/nb_iterations(2):(nb_iterations(1)-1/nb_iterations(2))];
%     semilogy(itNumbers,cat(2,err.ground_truth(:),err.l1(:),err.l2(:),err.l1Bregman(:),err.l2Bregman(:))); hold on;
%     semilogy(cat(1,err.ground_truth(end,:),err.l1(end,:),err.l2(end,:),err.l1Bregman(end,:),err.l2Bregman(end,:)).','*');
%     title('error'); xlabel('iteration'); legend({'\Delta ground truth','\Delta l1 target','\Delta l2 target','l1Bregman','l2Bregman'});
%     figure();
%     plot(observation, 'Color', [0.8, 0.0, 0.0], 'LineWidth', 2); hold on;
%     plot(x_optim, '-o', 'Color', [0, 0.6, 0], 'LineWidth', 2);
%     plot(ground_truth, 'Color', [0.1, 0.1, 0.1]);
%     legend({'Observation', 'splitBregman', 'Ground Truth'});
%
%     % find the peaks in 1D by minimizing the number of values that is non-zero
%     variableSize = [50, 1];  % A 1D problem
%     l1Weight = 100;  % how important is sparsity compared to best fit?
%     l1_subst_weight = 0.5;  % controls convergence rate
%     x0 = zeros(variableSize);
%     ground_truth = zeros(variableSize);
%     ground_truth([4]) = [1];
%     ground_truth = real(fft(ground_truth));  % A cosine
%     observation = ground_truth + 0.05*randn(variableSize);
%     l1_function = @(x) l1Weight*diff([x; 0]); % enforce sparsity of the gradient => blocks
%     l2_function = @(x) x-observation; % denoising, no deblurring
%     minimize_l2 = [];
%     nb_iterations = [50 10];  % The number of outer and inner iterations to do
%     stop_criterion = @(outerIt, l1_error, l2_error) l1_error<1e-6 && l2_error<1e-6;  % Alternative to the above
%     [x_optim, err] = splitBregmanMinimize(x0, l1_function, l2_function, minimize_l2, l1_subst_weight, nb_iterations, stop_criterion, true, ground_truth);
%     itNumbers = 1/nb_iterations(2)+[0:1/nb_iterations(2):(nb_iterations(1)-1/nb_iterations(2))];
%     semilogy(itNumbers,cat(2,err.ground_truth(:),err.l1(:),err.l2(:),err.l1Bregman(:),err.l2Bregman(:))); hold on;
%     semilogy(cat(1,err.ground_truth(end,:),err.l1(end,:),err.l2(end,:),err.l1Bregman(end,:),err.l2Bregman(end,:)).','*');
%     title('error'); xlabel('iteration'); legend({'\Delta ground truth','\Delta l1 target','\Delta l2 target','l1Bregman','l2Bregman'});
%     figure();
%     plot(observation, 'Color', [0.8, 0.0, 0.0], 'LineWidth', 2); hold on;
%     plot(x_optim, '-o', 'Color', [0, 0.6, 0], 'LineWidth', 2);
%     plot(ground_truth, 'Color', [0.1, 0.1, 0.1]);
%     legend({'Observation', 'splitBregman', 'Ground Truth'});
%
%     % Sparsity in Fourier space
%     variableSize = [50, 1];  % A 1D problem
%     l1Weight = 1;  % how important is sparsity compared to best fit?
%     l1_subst_weight = 0.5;  % controls convergence rate
%     x0 = zeros(variableSize);
%     ground_truth = zeros(variableSize);
%     ground_truth([4, 10]) = [1, -1/3];  % The first two harmonics of a square wave
%     ground_truth = real(fft(ground_truth));  % Force it to be real
%     plot(ground_truth)
%     observation = ground_truth + 0.05*randn(variableSize);
%     function shrinkable = l1_function(x)
%         x_ft = fft(x);  % sparsity in Fourier space
%         shrinkable = l1Weight * [real(x_ft); imag(x_ft)];
%         % The previous will minimize the Euclidian norm for (real, imag), i.e. abs(x_ft).
%     end
%     l2_function = @(x) x-observation; % denoising
%     minimize_l2 = [];
%     nb_iterations = [50 10];  % The number of outer and inner iterations to do
%     stop_criterion = [];  %@(outerIt, l1_error, l2_error) l1_error<1e-6 && l2_error<1e-6;  % Alternative to the above
%     [x_optim, err] = splitBregmanMinimize(x0, @l1_function, l2_function, minimize_l2, l1_subst_weight, nb_iterations, stop_criterion, true, ground_truth);
%     itNumbers = 1/nb_iterations(2)+[0:1/nb_iterations(2):(nb_iterations(1)-1/nb_iterations(2))];
%     semilogy(itNumbers,cat(2,err.ground_truth(:),err.l1(:),err.l2(:),err.l1Bregman(:),err.l2Bregman(:))); hold on;
%     semilogy(cat(1,err.ground_truth(end,:),err.l1(end,:),err.l2(end,:),err.l1Bregman(end,:),err.l2Bregman(end,:)).','*');
%     title('error'); xlabel('iteration'); legend({'\Delta ground truth','\Delta l1 target','\Delta l2 target','l1Bregman','l2Bregman'});
%     figure();
%     plot(observation, 'Color', [0.8, 0.0, 0.0], 'LineWidth', 2); hold on;
%     plot(x_optim, '-o', 'Color', [0, 0.6, 0], 'LineWidth', 2);
%     plot(ground_truth, 'Color', [0.1, 0.1, 0.1]);
%     legend({'Observation', 'splitBregman', 'Ground Truth'});  
%    
%     % Denoise 3D by minimizing the gradients (diff3D)
%      variableSize = [5 5 5];
%     l1Weight = 10;
%     l1_subst_weight = 1/10;
%     x0 = zeros(variableSize);
%     ground_truth = zeros(variableSize); ground_truth(2:4,2:4,1:3) = 1;
%     observation = ground_truth+0.10*randn(variableSize);
%     diff3D = @(x) l1Weight*cat(4,diff(x([end,1:end],:,:,1),1,1),diff(x(:,[end,1:end],:,2),1,2),diff(x(:,:,[end,1:end],3),1,3));
%     l1_function = @(x) Shrinkable(diff3D(repmat(x,[1 1 1 3])),3); % sparsity of gradient
%     l2_function = @(x) x-observation; % denoising
%     minimize_l2 = [];
%     splitBregmanMinimize(x0, l1_function,l2_function, minimize_l2, l1_subst_weight, [50 1], [], true, ground_truth);
% 
%     % Deblur and denoise 3D by minimizing the gradients (diff3D)
%     variableSize = [5 5 5];
%     l1Weight = 10;
%     l1_subst_weight = 1/10;
%     x0 = zeros(variableSize);
%     ground_truth = zeros(variableSize); ground_truth(2:3, 2:4, 1:3) = 1;
%     h = @(x) 0.5*x+0.25*circshift(x, [1 0 0])+0.25*circshift(x, -[1 0 0]);
%     observation = h(ground_truth)+0.01*randn(variableSize);
%     diff3D = @(x) l1Weight*cat(4, diff(x([end, 1:end], :, :, 1), 1, 1), diff(x(:, [end, 1:end], :, 2), 1, 2), diff(x(:, :, [end, 1:end], 3), 1, 3));
%     l1_function = @(x) Shrinkable(diff3D(repmat(x, [1 1 1 3])), 3); % sparsity of gradient
%     l2_function = @(x) h(x)-observation; % blur
%     minimize_l2 = [];
%     splitBregmanMinimize(x0, l1_function, l2_function, minimize_l2, l1_subst_weight, [50 1], [], true, ground_truth);
%
%     % Deblur and denoise 3D with user L2 solver
%     variableSize = [5 5 5];
%     l1Weight = 100;
%     l1_subst_weight = 1/10;
%     x0 = zeros(variableSize);
%     ground_truth = zeros(variableSize); ground_truth(2:3, 2:4, 1:3) = 1;
%     h = @(x) 0.5*x+0.25*circshift(x, [1 0 0])+0.25*circshift(x, -[1 0 0]);
%     observation = h(ground_truth)+0.00*randn(variableSize);
%     diff3D = @(x) l1Weight*cat(4, diff(x([end, 1:end], :, :, 1), 1, 1), diff(x(:, [end, 1:end], :, 2), 1, 2), diff(x(:, :, [end, 1:end], 3), 1, 3));
%     l1_function = @(x) Shrinkable(diff3D(repmat(x, [1 1 1 3])), 3); % sparsity of gradient
%     l2_function = @(x) h(x)-observation; % blur
%     l1_functions_on_conj = @(x) diff3D(x.components{1}([1, end:-1:2], [1, end:-1:2], [1, end:-1:2], :));
%     diracInput = zeros(variableSize); diracInput(1, 1, 1) = 1;
%     otf = fftn(h(diracInput));
%     diffMatrix = l1_function(diracInput); clear diracInput;
%     fft3 = @(x) fft(fft(fft(x, [], 1), [], 2), [], 3);
%     diffMatrixFftSqd = sum(abs(fft3(diffMatrix.components{1})).^2, 4); clear diffMatrix;
%     minimize_l2 = @(l1_target, l2_target_offset, l1wt) ifftn((sum(conj(fft3(l1_functions_on_conj(l1_target*l1wt))), 4)*l1wt + ...
%         conj(otf).*fft3(l2_target_offset+observation))...
%         ./max(diffMatrixFftSqd*(l1wt^2) + abs(otf).^2, eps(1)), 'symmetric');
%     splitBregmanMinimize(x0, l1_function, l2_function, minimize_l2, l1_subst_weight, [100 4], [], true, ground_truth);
%
function [x_optim, err] = splitBregmanMinimize(x0, l1_function, l2_function, minimize_l2, l1_substitute_convergence_weight, nb_iterations, stop_criterion, want_statistics, ground_truth, progress_callback)
    %
    % Initialize defaults for testing
    %
    if nargin<1 && nargout<1
        close all;
    end    
    if nargin<1
        x0 = zeros([1 1]*750);
    end
    % Regularization specific
    if nargin<2 || isempty(l1_function)
        if ~isnumeric(x0)
            error('l1_function should be specified if x0 is not a matrix.');
        end
            
        regularizationWeight = 10;
        
        diff_x = @(x) diff(x([end,1:end],:),1,1);
        diff_y = @(x) diff(x(:,[end,1:end]),1,2);
        l1_functions_parallel = @(x) regularizationWeight*cat(3,diff_x(x(:,:,1)),diff_y(x(:,:,2)));
        l1_function = @(x) Shrinkable(l1_functions_parallel(repmat(x,[1 1 2])),2);
        clear diff_x diff_y regularizationWeight;
        l1_functions_on_conj = @(x) l1_functions_parallel(x.components{1}([1,end:-1:2],[1,end:-1:2],:)); % also used by minimize_l2
        clear l1_functions_parallel;
        
        % diffMatrixFftSqd also required by the default for minimize_l2
        diracInput = x0; diracInput(1,1) = 1;
        diffMatrix = l1_function(diracInput);
        diffMatrixFftSqd = sum(abs(fft2(diffMatrix.components{1})).^2, 3); % Also required by the default for minimize_l2
        clear diracInput diffMatrix;
        
        %
        % Define the related input arguments now too
        %
        [observation, ground_truth, otf, noise_level] = getExampleObservation(size(x0));
        l2_function = @(x) ifft2(otf.*fft2(x),'symmetric')-observation;
        % ground_truth is now know, so it can be used for diagnostics
        minimize_l2 = @(l1_target,l2_target_offset,l1wt) ifft2((sum(conj(fft2(l1_functions_on_conj(l1_target*l1wt))),3)*l1wt + ...
                conj(otf).*fft2(l2_target_offset+observation))...
                ./max(diffMatrixFftSqd*(l1wt^2) + abs(otf).^2,eps(1)),'symmetric');
        clear noise_level diffMatrixFftSqd otf observation;
    end
    if ~exist('l2_function','var') || isempty(l2_function)
        error('If the l1_function arguments is defined, so must the l2_function argument.');
    end
    if nargin<5 || isempty(l1_substitute_convergence_weight)
        l1_substitute_convergence_weight = 0.10;  % Make sure that l1-norm converges before the l1_substitute
    end
    if nargin<6 || isempty(nb_iterations)
        nb_iterations = [30 3];
    end
    nb_iterations(end+1:3) = 1;
    if nargin<7 || isempty(stop_criterion)
        stop_criterion=@(outer_iter,l1_error,l2_error) false;
    end
    if nargin<8 && ~exist('want_statistics', 'var')
        want_statistics = (nargout>1);
    end
    if nargin<9 && ~exist('ground_truth', 'var')
        ground_truth = [];
    end
    if nargin<10 || isempty(progress_callback)
        progress_callback = [];
    end
    
    if ~exist('minimize_l2','var') || isempty(minimize_l2)
        minimize_l2=@(l1_target, l2_target_offset, l1_subst_weight) ...
            genericMinimizeL2(size(x0), nb_iterations(3),...
                @(x) l1_subst_weight*(l1_function(x)-l1_target),... % l1_subst_weight under the l2-norm
                @(x) l2_function(x)-l2_target_offset);
    end
    
    %
    % Do optimization
    %
    % Initialize l1 substitute and Bregman constrainers
    l1_substitute = l1_function(x0)*0; % Get the size(s) of the L1 vector(s)
    l1_Bregman_constrainers = l1_substitute;
    l2_Bregman_constrainer = l2_function(x0)*0; % Get the size of the L2 vector
    x_optim = x0;
    clear x0
    
    % Initialize an error struct for analysis
    if ~isempty(ground_truth)
        err  =struct();
        err.l1 = [];
        err.l2 = [];
        err.l1Bregman = [];
        err.l2Bregman = [];
        err.ground_truth = [];
    end
    % Constrain the unconstrained problem
    for outer_iter = 1:nb_iterations(1)
        % Solve unconstrained problem:
        % min_u sum(abs(l1_function(u))) + sum(abs(l2_function(u)-(0-l2_Bregman_constrainer)).^2)
        for inner_iter = 1:nb_iterations(2)
            % Minimize the l2-norms of l2_function and the distance l1_function-l1_substitute:
            % min_u l1_substitute_convergence_weight^2*sum(abs(l1_function(u)-(l1_substitute-l1_Bregman_constrainers)).^2) + ...
            %                                     sum(abs(l2_function(u)-(0-l2_Bregman_constrainer)).^2);
            all_min_l2_args = {l1_substitute-l1_Bregman_constrainers, 0-l2_Bregman_constrainer, l1_substitute_convergence_weight, x_optim, nb_iterations(3)};
            all_min_l2_args = {all_min_l2_args{1:min(end,nargin(minimize_l2))}};
            x_optim = minimize_l2(all_min_l2_args{:});
            clear all_min_l2_args;
            
            % The weight of the distance l1_function-l1_substitute is increased first by the l1BregmanConstrainer
            
            % Minimize the l2-distance of d to l1_substitute while keeping the l1-norm of d minimized:
            % min_d sum(abs(d)) + l1_substitute_convergence_weight^2*sum(abs(l1_function_current_vectors+l1_Bregman_constrainers - d).^2)
            l1_function_current_vectors = l1_function(x_optim);
            l1_substitute = shrink(l1_function_current_vectors+l1_Bregman_constrainers,1.0/l1_substitute_convergence_weight^2);

            % update Bregman parameter for L1 optimization
            l1_Bregman_constrainers = l1_function_current_vectors + (l1_Bregman_constrainers-l1_substitute);
            
            if ~isempty(progress_callback)
                try
                    all_input_args = {[outer_iter, inner_iter], nb_iterations, x_optim, l1_substitute, l1_Bregman_constrainers, l2_Bregman_constrainer};
                    all_input_args = {all_input_args{1:min(end,nargin(progress_callback))}};
                    progress_callback(all_input_args{:});
                    clear all_input_args
                catch Exc
                    fprintf('Progress report callback function progress_callback failed!\n\nError message details:\n%s\n',Exc.getReport('extended'));
                end
            end
            
            l1_error = mean(abs(l1_function_current_vectors(:)));
            clear l1_function_current_vectors;
            
            % Store errors in error struct if requested
            if want_statistics
                if ~isempty(ground_truth)
                    if ~isstruct(ground_truth)
                        err.ground_truth(inner_iter,outer_iter) = sqrt(mean(abs(x_optim(:)-ground_truth(:)).^2));
                    else
                        % In case the ground truth is in the form of a struct
                        fieldNames = sort(fields(ground_truth));
                        ground_truth_error = 0;
                        for fieldIdx = 1:numel(fieldNames)
                            ground_truth_field = getfield(ground_truth,fieldNames{fieldIdx});
                            x_optim_field = getfield(x_optim,fieldNames{fieldIdx});
                            ground_truth_error = ground_truth_error+sqrt(mean(abs(x_optim_field(:)-ground_truth_field(:)).^2));
                        end
                        err.ground_truth(inner_iter,outer_iter) = ground_truth_error;
                        clear fieldNames x_optim_field ground_truth_field ground_truth_error
                    end
                else
                    err.ground_truth(inner_iter,outer_iter) = Inf;
                end
                err.l1(inner_iter,outer_iter) = l1_error;
                err.l1Bregman(inner_iter,outer_iter) = sqrt(mean(abs(l1_Bregman_constrainers(:)).^2));
                err.l2Bregman(inner_iter,outer_iter) = sqrt(mean(abs(l2_Bregman_constrainer(:)).^2));
            end
        end
        
        % Add errorL2Target back to the observation to iteratively solve the constrained problem
        l2_function_current_vector = l2_function(x_optim);
        l2_Bregman_constrainer = l2_function_current_vector + l2_Bregman_constrainer;
        
        l2_error = sqrt(mean(abs(l2_function_current_vector(:)).^2));
        clear l2_function_current_vector;
        
        if want_statistics
            for inner_iter = 1:nb_iterations(2)
                err.l2(inner_iter,outer_iter) = l2_error;
            end
        end
        
        % Check if the stop criterion is reached
        if outer_iter<nb_iterations(1)
            input_args = {outer_iter,l1_error,l2_error};
            if stop_criterion(input_args{1:nargin(stop_criterion)})
                break;
            end
        end
    end
    clear l1_Bregman_constrainers l1_substitute l2_Bregman_constrainer;
    
    %
    % Output
    %
    if want_statistics && nargout<2
        figure();
        nb_outer_it = size(err.l2,2);
        itNumbers = 1/nb_iterations(2)+[0:1/nb_iterations(2):(nb_outer_it-1/nb_iterations(2))];
        semilogy(itNumbers, cat(2,err.ground_truth(:),err.l1(:),err.l2(:),err.l1Bregman(:),err.l2Bregman(:))); hold on;
        semilogy(cat(1,cat(1,err.ground_truth(end,:),err.l1(end,:),err.l2(end,:),err.l1Bregman(end,:),err.l2Bregman(end,:)).',NaN*ones(1,5)),'*');
        title('error');
        xlabel('iteration');
        legend({'\Delta ground truth','\Delta l1 target','\Delta l2 target','l1Bregman','l2Bregman'});
        fprintf('Difference with ground truth: %0.8f%%\n',100*err.ground_truth(end));
    end
    if nargout<1 && isnumeric(x_optim)
        figure();
        switch(ndims(x_optim))
            case {1,2}
                if any(size(x_optim)==1)
                    plot(x_optim);
                    xlabel('x'); ylabel('intensity');
                else
                    imagesc(real(x_optim));
                    xlabel('x'); ylabel('y');
                    axis image equal tight;
                end
            otherwise
                subplot(2,2,1);
                imagesc(x_optim(:,:,1,1).');
                xlabel('x'); ylabel('y'); title('z: front'); axis image equal tight;
                colorbar();
                subplot(2,2,2);
                imagesc(x_optim(:,:,1+floor(end/2),1).'); axis image equal tight;
                xlabel('x'); ylabel('y'); title('z: center');
                colorbar();
                subplot(2,2,3);
                imagesc(real(squeeze(x_optim(:,1,:,1)).')); axis image equal tight;
                xlabel('x'); ylabel('z'); title('y: top');
                colorbar();
                subplot(2,2,4);
                imagesc(real(squeeze(x_optim(:,1+floor(end/2),:,1)).')); axis image equal tight;
                xlabel('x'); ylabel('z'); title('y: middle');
                colorbar();
        end
        clear x_optim;
    end
end


%
% Auxiliary functions in case not all input arguments are given
%

% Get an example image of a given size (padded 750x750px image)
function [observation, ground_truth, otf, noise_level] = getExampleObservation(img_shape)
    noise_level = 0.01;

%     ground_truth = 1-getTestImage('usaf1951_1500x1500.png');
    ground_truth = 1-getTestImage('usaf');
    img_shape(end+1:3) = 1;
    ground_truth(1+img_shape(1),1+img_shape(2),1+img_shape(3)) = 0;
    ground_truth = ground_truth(1:img_shape(1),1:img_shape(2),1:img_shape(3));
    
    randn('seed', 0);
    noise = noise_level*randn(size(ground_truth));

    psf = fspecial('gaussian', size(ground_truth), 5);
    otf = fft2(ifftshift(psf));
    clear psf;
    % otf is also required by the default for l2_function and minimize_l2
    observation = ifft2(otf.*fft2(ground_truth),'symmetric')+noise;
    clear noise;
end

% Optimizes u for sum_i ||fi(u)||_2^2
%
% Input:
%     varSize: the matrix size of the variable u to optimize
%     nb_iterations: number of iterations to perform
%     fi(u): one or more functions to be squared and summed
%
function optimX = genericMinimizeL2(varargin)
    varSize = varargin{1};
    if prod(varSize) > 1000
      logMessage('Number of variables is large (%d), provide the minimize_l2 function to make this efficient.', prod(varSize));
    end
    nb_iterations = varargin{2};
    functionsToMinimize = varargin(3:end);
    function result = optimFunc(x)
        result = [];
        for argIdx = 1:numel(functionsToMinimize)
            term = functionsToMinimize{argIdx}(reshape(x, varSize));
            result(end+[1:numel(term)]) = term(:);
        end
    end

    optimX = lsqnonlin(@optimFunc, zeros(prod(varSize), 1), [], [], optimset('MaxIter', nb_iterations, 'Display', 'off'));
    
    optimX = reshape(optimX,varSize);
end
