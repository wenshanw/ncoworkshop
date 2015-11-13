for yyyy in {2008..2013}; do
  echo ${yyyy}

  cp map2stn.ncl map2stn_${yyyy}.ncl

  sed -i "s/yr=\"2010\"/yr=\"${yyyy}\"/g" ${exe}_${yyyy}.ncl

  ncl ${exe}_${yyyy}.ncl&>out_${yyyy}&
done
