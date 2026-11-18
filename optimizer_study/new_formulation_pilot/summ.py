import json,sys,glob,pandas as pd
pd.set_option("display.width",220,"display.max_rows",400,"display.max_columns",30)
r=pd.concat([pd.read_json(f,lines=True) for f in glob.glob(sys.argv[1])])
g=r.groupby(["problem","arm","optimizer"]).agg(n=("E_end","size"),E_start=("E_start","mean"),E_end=("E_end","mean"),E_min=("E_end","min"),ground=("E_ground","first"),fb_start=("P_fb_start","mean"),fb_end=("P_fb_end","mean"),fb_max=("P_fb_end","max"),raw_end=("P_raw_end","mean"),raw_max=("P_raw_end","max"),evals=("evals","median")).round(3).reset_index()
print(g.sort_values(["problem","arm","E_end"]).to_string(index=False))
