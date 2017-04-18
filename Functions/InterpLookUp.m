function [ InterpTable ] = InterpLookUp( number_slices, aspect_ratio )

% INTERPLOOKUP Small function that finds the mapping of the slice number to the
% corresponding z-coordinate in the interpolated volume.

%Inputs: 
    % number_slices: number of slices in the given sample
    % aspect_ratio: z-step / x-y step
    
% Output:
    % InterpTable: Look-up table in the form of a matrix. The first column
    % holds the slice number and the 2nd column holds the new z-coordinate.

% 10.07.2016 A.Harel

M=zeros(number_slices,number_slices);
for i=1:number_slices
    for j=1:number_slices
        if i==j
            M(i,j)=1;
        else
        end
    end
end

M=double(M);
ny=size(M,2);nx=size(M,1)*aspect_ratio; %% desired output dimensions
[xq_test,yq_test]=  ndgrid(linspace(1,size(M,1),nx),...
    linspace(1,size(M,2),ny));
Mint=interp2(M,xq_test,yq_test);
Mint=Mint>=0.9;

[InterpTable(:,2),InterpTable(:,1)]= ind2sub(size(Mint),find(Mint));

end
