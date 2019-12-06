
for ii = -Ng:Ng
    for j = -Ng:Ng
        r2 = sqrt(ii^2 + j^2);
        if r2 <= Ng
          if r2 ~= 0
              if ii < 0 
                 ii = ii + 2*Ng + 1;
              end
              if j < 0 
                 j = j + 2*Ng + 1;
              end
              ii = ii + 1;
              j = j + 1;
              Vfft(ii,j) = 20/(norm([ii,j],2)^2 + 1);
          else
              Vfft(ii,j) = 10;
          end  
        end
    end
end