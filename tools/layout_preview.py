"""Render actual exported example GDS polygons (requires [layout,plot])."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch,Patch
import klayout.db as db

from twpa.layout import load_process
from twpa.layout.generate import generate,KITLayout,JosephsonLayout


def render(output):
    colors={'M1':'#196a9d','M2':'#e49b3e','VIA':'#29984a','JJ':'#db384e','CAP':'#9a67bc','IND':'#1aa99a'}
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    fig.set_facecolor('#fafbfc')
    families=['kit','jpa','squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm']
    titles=['Research KIT · two-row excerpt','Resonant JPA · illustrative','SQUID JPA · illustrative',
            '4WM JTWPA · illustrative','rf-SQUID 3WM · illustrative','SNAIL 3WM · illustrative']
    folder=Path('validation/runs/layout/previews'); folder.mkdir(parents=True,exist_ok=True)
    for family,title,ax in zip(families,titles,axes.flat):
        p=load_process('kit_nycreates_research' if family=='kit' else 'jj_overlap_illustrative')
        d=KITLayout(rows=2,periods_per_row=2) if family=='kit' else JosephsonLayout(family=family,cells=1 if 'jpa' in family else 2)
        path=folder/(family+'.gds'); generate(d,p,path)
        ly=db.Layout(); ly.read(str(path)); top=ly.top_cell()
        roles=['M1'] if family=='kit' else list(colors)
        bb=None
        for role in roles:
            reg=db.Region(top.begin_shapes_rec(ly.layer(*p.layer(role)))).merged()
            if not reg.is_empty(): bb=reg.bbox() if bb is None else bb+reg.bbox()
            for poly in reg.each():
                vertices=[]; codes=[]
                for ring in [list(poly.each_point_hull())]+[list(poly.each_point_hole(i)) for i in range(poly.holes())]:
                    vertices.extend([(v.x*ly.dbu,v.y*ly.dbu) for v in ring]+[(ring[0].x*ly.dbu,ring[0].y*ly.dbu)])
                    codes.extend([MPath.MOVETO]+[MPath.LINETO]*(len(ring)-1)+[MPath.CLOSEPOLY])
                ax.add_patch(PathPatch(MPath(vertices,codes),facecolor=colors[role],edgecolor='none',alpha=.9))
        b=bb.to_dtype(ly.dbu); pad=max(b.width(),b.height())*.05
        ax.set(xlim=(b.left-pad,b.right+pad),ylim=(b.bottom-pad,b.top+pad),xlabel='x (um)',ylabel='y (um)',title=title)
        ax.set_aspect('equal',adjustable='box'); ax.tick_params(labelsize=8)
        ax.spines[['top','right']].set_visible(False)
    fig.legend(handles=[Patch(color=color,label=role) for role,color in colors.items()],loc='outside lower center',ncol=6,frameon=False)
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=150); plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='docs/assets/layout-families.png')
    render(p.parse_args().output)
