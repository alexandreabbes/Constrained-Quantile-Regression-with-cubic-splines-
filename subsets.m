%here is a function finding all subsets of a set
%subsets
function subA=subsets(A)
 
 bitNo = length(A);     % number of bits equals to the number of input elements
 setNo = 2 ^ bitNo ; % number of sets 
 subA=cell(setNo,1);
 if bitNo==0
     subA={[]}
 else
     
  for setId = 1 : setNo
   % convert number to a binary string and that to logical indices
   setIndices = 1+bitNo-find(dec2bin(setId,bitNo)==dec2bin(1));
   % select the current set by using the logical indices
   currentSet = A(setIndices);
  
   subA{setId+1}=[currentSet];
   %disp(currentSet)
  end
 end
 
 end
