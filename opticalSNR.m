function oSNR = opticalSNR( d )
v_d=std(d);
m_d=mean(d);
oSNR=m_d./v_d;