createdb dsi
psql -c 'create extension rdkit' dsi
psql -c 'create table tot_data (id SERIAL, smiles text, enamine_id text,parent_id text)' dsi
cat DSI_new.smi | sed '1d; s/\\/\\\\/g' | psql -c "copy tot_data (smiles,enamine_id,parent_id) from stdin with delimiter ' '" dsi

psql dsi - then
select * into mols from (select id,enamine_id,mol_from_smiles(smiles::cstring) m from tot_data) tmp where m is not null;
create index molidx on mols using gist(m);
select enamine_id,torsionbv_fp(m) as torsionbv,morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into fps from mols;
create index fps_ttbv_idx on fps using gist(torsionbv);
create index fps_mfp2_idx on fps using gist(mfp2);
create index fps_ffp2_idx on fps using gist(ffp2);
create or replace function get_mfp2_neighbors(smiles text)
    returns table(id integer, m mol, similarity double precision) as
  $$
  select enamine_id,m,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
  from fps join mols using (enamine_id)
  where morganbv_fp(mol_from_smiles($1::cstring))%mfp2
  order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
  $$ language sql stable;