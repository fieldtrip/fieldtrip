function plotCov(N,s)


F = zeros(length(s),N);

for j=1:length(s)
	K = construct_prior(N,-s(j));
	K = scale_prior_bc(K,'lambda',1);
	V = inv(full(K));
	F(j,:) = V(1,:);
end

figure;
hold on;
	for j=1:length(s)
		plot(1:N,F(j,:));
	end
hold off;