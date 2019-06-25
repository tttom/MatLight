%
% password = createPassword(len)
%
function password = createPassword(len)
  rng = RandStream('mt19937ar', 'Seed', mod(floor(now() * 1e8), 2^32));
  
  character_str = ['A':'Z', 'a':'z', '0':'9', ':.,<>?/[]-_+=#@;*!$%&'];
  
  password = character_str(rng.randi(numel(character_str), 1, len));
end