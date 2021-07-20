function B_mat=bs_function(degree,interior_knots,boundary_knots)

boundary_knots=sort(boundary_knots);
interior_knots_sorted=[];
if ~isempty(interior_knots)
    interior_knots_sorted=interior_knots;
end
knots=[repelem(boundary_knots(1),degree+1),interior_knots_sorted,repelem(boundary_knots(2),degree+1)];
K=length(interior_knots)+degree+1;
for j=1:K
    B_mat{j}=basis_function(degree,j,knots,K);
end


% for i=1:length(x)
%     if x(i)==boundary_knots(2)
%         B_mat(i,K)=1;
%     end
% end
% if intercept==0
%     B_mat=B_mat(:,1:K-1);
% end
    