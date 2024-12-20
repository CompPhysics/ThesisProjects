% Constants for self-consistent solving of in-medium self-energy Sigma and
% vertex function Gamma using Green's functions formalism

% hbar times c in units of MeV
hbarc = 197.327053;

% Proton and neutron masses and reduced mass, in units of MeV
m_p = 938.27231;
m_n = 939.56563;
redmass_np = (m_p.*m_n)./(m_p+m_n);

% hbar*c square divided by the proton mass
meq1 = (hbarc.*hbarc)./m_p;
