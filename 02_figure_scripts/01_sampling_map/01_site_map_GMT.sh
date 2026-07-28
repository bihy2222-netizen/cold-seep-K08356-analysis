#!/bin/bash
# 恢复的 site_map_2d 生成脚本

cd "/Users/catherine/Downloads/GMT_demo"
mkdir -p output
cd output

echo "=== 恢复的 site_map_2d 生成脚本 ==="

# 数据文件路径
DATA_FILE="../gebco_2024/GEBCO_2024.nc"

if [ ! -f "$DATA_FILE" ]; then
    echo "错误: 找不到数据文件"
    exit 1
fi

echo "使用数据: $(basename $DATA_FILE)"

# 创建配色
gmt makecpt -Crainbow -T0/3000/100 > map.cpt

echo "生成 site_map_2d.pdf..."

gmt begin site_map_2d pdf
    # 基本设置
    gmt set MAP_FRAME_TYPE fancy
    
    # 区域设置
    RANGE="110/112/16/18"
    
    # 绘制地形
    gmt grdimage $DATA_FILE -R$RANGE -JM15c -Cmap.cpt
    gmt coast -W0.5p,black -Da
    
    # 绘制站点
    gmt plot -Sc0.25c -Gred -W0.5p,white << EOF
111.11707 17.718045
111.12198 17.7036
111.12114 17.70335
111.114929 17.699036
111.12076 17.70335
111.05574 17.62273
110.4 16.9
110.47 16.73
110.4 16.69
110.4724 16.7283
110.4719 16.7285
110.4718 16.7284
110.43 16.71
110.5 16.5
111 17.2
111 17
EOF
    
    # 地图装饰
    gmt basemap -Baf -BWSen+t"Site Map 2D"
    gmt colorbar -Cmap.cpt -DjBC+w10c/0.4c -Bx+l"Elevation (m)"
    gmt basemap -TdjTR+w1.5c+f
    
gmt end show

rm -f map.cpt

echo "✓ site_map_2d.pdf 已生成"
