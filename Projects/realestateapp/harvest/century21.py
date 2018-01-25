import requests
from bs4 import BeautifulSoup
import pandas

# FOR WEB CRAWLING, MAKE FOR LOOP WITH REQUEST AND STORE CSV
r = requests.get("https://www.century21.com/real-estate/rock-springs-wy/LCWYROCKSPRINGS/")
content = r.content
soup = BeautifulSoup(content, "html.parser")
harvest = soup.find_all("div", {"class": "property-card-primary-info"})
l = []
for index in harvest:

#list_price = harvest[0].find("a").text.replace("\n", "").replace(" ", "")
    d = {}
    list_price = index.find("a").text.replace("\n", "").replace(" ", "")
    try:
        list_price = index.find("a").text.replace("\n", "").replace(" ", "")
        beds = index.find("div", {"class": "property-beds"})
        beds = beds.find("strong").text
    
        sqft = index.find("div", {"class": "property-sqft"})
        sqft = sqft.find("strong").text

        d["list_price"] = list_price 
        d["beds"] = beds
        d["sqft"] = sqft
         
    except:
        pass
    l.append(d)

    print(list_price, beds, sqft)
print(l)
df = pandas.DataFrame(l)
print(df)
df.to_csv("century.csv")



