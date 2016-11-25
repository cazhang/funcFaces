 function Vout = mapClosestDeltaMutu(X, Vin, Nsets, n)
 
 fprintf('Re-ordering basis: shape 1 as reference. \n');
 k = -n;
 
 T = M2T(X);
 Vout = circshift(Vin, k, 3);
 for i=2:Nsets
     V1 = Vout(:,:,1);
     V2 = Vout(:,:,i);
     %M = X(:,:,i-1)';
     idx1 = mod(n+1, Nsets);
     idx2 = mod(n+i, Nsets);
     if idx1 == 0
         idx1 = Nsets;
     end
     if idx2 == 0
         idx2 = Nsets;
     end
     M = T(:, :, idx2, idx1);
     R = closestRotation(M);
     Vc =  R * V1;
     dest = annquery(V2, Vc, 1);
     
     Vout(:,:,i) = V2(:, dest);
 end
 
        
 end