function Conc_CH4 = Conc_CH4_Model_DOM_Aug (p, x)

Conc = Along_Canal_Model_DOM_Aug (p,x);
Conc_CH4 = Conc(:,3)';

end

