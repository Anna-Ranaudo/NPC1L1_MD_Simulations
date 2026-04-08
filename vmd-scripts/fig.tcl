display resize 3000 2000
display projection Orthographic

# luci migliori
light 0 on
light 1 on
light 2 off
light 3 off

display shadows on
display ambientocclusion on
display aoambient 0.85
display aodirect 0.35


render TachyonInternal figura.tga

ffmpeg -i figura.tga figura.png
