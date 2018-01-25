import folium
import pandas as pd

SF_CORD = [37.76, -122.45]
MAX_RECORDS = 1000

crime_data = pd.read_csv('./Data/Crime.csv')

crime_map = folium.Map(location=SF_CORD, zoom_start=12)
crime_map.save("SFcrime.html")

for each in crime_data[0:MAX_RECORDS].iterrows():
    folium.Marker(location = [each[1]['Y'],each[1]['X']]).add_to(crime_map)


crime_map.save("SFcrime.html")
