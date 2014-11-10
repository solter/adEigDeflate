#makes a spiked matrix with spikes
function [S] = makeS(n, embedded, nz = n-2, randPos = true)

  S = schur(randn(n-2));
  S = [randn(n-2,1), S, randn(n-2,1)];
  S = [randn(1,n); S; 0, randn(1,n-1)];

  kt = n-1; kb = 2;
  for i=1:nz #make m zeros (the theoretical minimum number of spike tips to deflate
    if(randPos == true)
      zeroed = false;
      while(~zeroed)
        
        k = 1 + ceil((n-2)*rand(1));
        if(rand(1) > .5)
          if(S(k,1) ~= 0)
            S(k,1) = 0;
            zeroed = true;
          end
        else
          if(S(end,k) ~= 0)
            S(end,k) = 0;
            zeroed = true;
          end
        end
      end
    else
      if(n - kt < kb )
        S(kt--,1) = 0;
      else
        S(end,kb++) = 0;
      end
    end
  end

  if(embedded)
    S2 = hessen(randn(3*n));
    S2(n+1:2*n,n+1:2*n) = S;
    S = S2;
  end

end
