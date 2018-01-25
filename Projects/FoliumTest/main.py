import folium

map = folium.Map(location=[38.58, -99.09], zoom_start=6)
fg = folium.FeatureGroup('Markers')

mark1 = folium.Marker(location=[38.59, -99.09], popup="HELLO",
            icon=folium.Icon(color='green'))

fg.add_child(mark1)

map.add_child(fg)
map.save("Amap.html")
