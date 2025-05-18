load {{ fn }} FILTER "nosort"
set showhydrogens True;

{{ orient }}

function _setModelState() {

  select;
  {% if not args.nospacefill %}Spacefill 0.0;{% endif %}

  frank off;
  font frank 16.0 SansSerif Plain;
  select *;
  set fontScaling false;

}

_setModelState;

background white
mo fill
mo translucent
mo cutoff {{ config.Cutoff }}

mo nomesh
mo COLOR {{ config.ColorNegative }} {{ config.ColorPositive }}
mo resolution {{ config.Resolution}}
mo titleformat " "

{% for mo, mo_fn in mos_fns %}
mo {{ mo }}
write image pngt "{{ mo_fn }}"
print "Wrote {{ mo_fn }}"
{% endfor %}
