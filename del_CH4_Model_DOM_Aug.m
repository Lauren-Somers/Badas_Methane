function del_CH4 = del_CH4_Model_DOM_Aug (p, x)

Conc = Along_Canal_Model_DOM_Aug (p,x);
del_CH4 = Conc(:,4)';

end
