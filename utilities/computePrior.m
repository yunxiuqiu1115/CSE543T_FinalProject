function a = computePrior(pvi, experienced, republican)
    democratic = 1 - republican;
    a = 18.98808 + pvi*0.66579 + experienced*5.63548 + republican*23.98806 + democratic*25.74237 - pvi*republican*1.2485;
    a = a/100;
end