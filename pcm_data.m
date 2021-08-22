

phi = [0.1];
rho_pcm = [770];
a = rho_npcm(phi,rho_pcm)

function x = rho_npcm(phi,rho_pcm)
    rho_np = 3600;
  x = rho_np.*phi + rho_pcm.*(1-phi);                    % function calculates the value of pho_npcm
end
