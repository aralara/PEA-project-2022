% Ferro Lara 10622035 - QR code ID 78143046 - Project C
clear all;
clc;

% Fitting distributions

trace = ["TraceC-Spec.txt" "TraceC-Design.txt" "TraceC-Breadbrd.txt" "TraceC-Software.txt" "TraceC-Test.txt"];

for trace_idx=1:5 % for each trace
	path = "./Traces/"+trace(1,trace_idx);
	Ai = importdata(path); 
	AiSorted = sort(Ai);
	N = size(Ai(:,1),1);

	if (trace_idx <= 2)
		t = [0:8000]/10;
	elseif (trace_idx == 3)
		t = [0:35000]/10;
	elseif (trace_idx == 4)
		t = [0:60000]/10;
	else
		t = [0:14000]/10;
	end

	% options for fsolve function
	options = optimoptions('fsolve');  
	options.MaxIterations = 1000; 
	options.MaxFunctionEvaluations = 5000;
	options = optimset('Display','off');

	% METHOD OF MOMENTS
	fprintf(1, "METHOD OF MOMENTS %s\n", trace(1,trace_idx));

	m1 = sum(Ai(:, 1))/N; 
	m2 = sum(Ai(:, 1) .^2)/N;
	m3 = sum(Ai(:, 1) .^3)/N;
	cv = (sqrt(sum((Ai(:,1)-m1) .^2)/N))/m1;
	fprintf(1, "First moment: %g\n", m1);
	fprintf(1, "Second moment: %g\n", m2);
	fprintf(1, "Third moment: %g\n", m3);
	fprintf(1, "Coefficient of variation: %g\n\n", cv);

	%	Uniform 
	a_unif = m1 - 0.5*sqrt(12*(m2 - m1.^2));
	b_unif = m1 + 0.5*sqrt(12*(m2 - m1.^2));
	y_unif = max(0, min(1, (t - a_unif)/(b_unif - a_unif)));
	fprintf(1, "Uniform left bound: %g\n", a_unif);
	fprintf(1, "Uniform right bound: %g\n\n", b_unif);

	%	Exponential
	lambda_exp = 1 / m1;
	y_exp = 1 - exp(-lambda_exp*t);
	fprintf(1, "Exponential rate: %g\n\n", lambda_exp);

	%	Two stages hyper-exponential (not pox cv < 1)
	if (cv > 1)
		params = [0.8/m1, 1.2/m1, 0.4];
		HyperExpMoments = @(p) [p(1,3)/p(1,1)+(1-p(1,3))/p(1,2), 2*(p(1,3)/p(1,1)^2+(1-p(1,3))/p(1,2)^2), 6*(p(1,3)/p(1,1)^3+(1-p(1,3))/p(1,2)^3)];
		F = @(p) (HyperExpMoments(p) ./ [m1, m2, m3]-1);
		res = fsolve(F, params, options);
		HyperExpCDF = @(p, t) 1-(p(1,3)*exp(-p(1,1)*t))-((1-p(1,3))*exp(-p(1,2)*t));
		y_hyperExp = HyperExpCDF(res, t);
		fprintf(1, "Hyper-Exponential rate one: %g\n", res(1,1));
		fprintf(1, "Hyper-Exponential rate two: %g\n", res(1,2));
		fprintf(1, "Hyper-Exponential probability one: %g\n\n", res(1,3));
	else 
		fprintf(1, "No solution for Hyper-Exponential: cv < 1\n\n");
	end

	%	Two stage hypo-exponential (not pox cv > 1)
	if (cv < 1)
		params = [1/(0.3*m1), 1/(0.7*m1)];
		HypoExpMoments = @(p) [1/p(1,1)+1/p(1,2), 2*(1/(p(1,1)^2)+1/(p(1,1)*p(1,2))+1/(p(1,2)^2))];
		F = @(p) (HypoExpMoments(p) ./ [m1, m2]-1);
		res = fsolve(F, params, options);
		HypoExpCDF = @(p, t) 1-(p(1,2)*exp(-p(1,1)*t)/(p(1,2)-p(1,1)))+(p(1,1)*exp(-p(1,2)*t)/(p(1,2)-p(1,1)));
		y_hypoExp = HypoExpCDF(res, t);
		fprintf(1, "Hypo-Exponential rate one: %g\n", res(1,1));
		fprintf(1, "Hypo-Exponential rate two: %g\n\n", res(1,2));
	else 
		fprintf(1, "No solution for Hypo-Exponential: cv > 1\n\n");
	end
	
	%	Erlang
	if (cv < 1)
		k_erlang = round(1/(cv^2));
		lambda_erlang = k_erlang / m1;
		
		sum_erlang = 0;
		for n_erlang = 0:(k_erlang-1)
			sum_erlang = sum_erlang + ((1/factorial(n_erlang))*exp(-lambda_erlang*t).*(lambda_erlang*t).^n_erlang);
		end
		y_erlang = 1 - sum_erlang;
		
		fprintf(1, "K: %g\n", k_erlang);
		fprintf(1, "Erlang rate: %g\n\n", lambda_erlang);
	end
	
	%	Compare the CDF of the distributions with the CDF obtained directly from the dataset
	f1 = figure('Name',trace(1,trace_idx));
	if (cv > 1) % Ok with hyper-exponential
		p = plot(AiSorted(:,1), [1:N]/N, "+", t, y_unif, "-", t, y_hyperExp, "-", t, y_exp, "-");
		legend("Dataset", "Uniform", "Two stage hyper-exponential", "Exponential");
	else %(cv < 1) Ok with hypo-exponential, erlang
		p = plot(AiSorted(:,1), [1:N]/N, "+", t, y_unif, "-", t, y_hypoExp, "-", t, y_exp, "-", t, y_erlang, "-");
		legend("Dataset", "Uniform", "Two stage hypo-exponential", "Exponential", "Erlang"); 
		p(5).LineWidth = 2;
	end
	p(2).LineWidth = 2; p(3).LineWidth = 2; p(4).LineWidth = 2;
	
	fprintf(1, "\n\n");
end