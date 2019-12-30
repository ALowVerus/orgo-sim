import requests
from bs4 import BeautifulSoup as bs
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By


def get_rendered(url, wait_key):
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('headless')
    w = webdriver.Chrome(options=chrome_options)
    w.get(url)
    delay = 3  # seconds
    my_elem = WebDriverWait(w, delay).until(EC.presence_of_element_located(wait_key))
    html = w.find_element(By.TAG_NAME, "html").get_attribute("innerHTML")
    w.quit()
    return bs(html, "html.parser")


def get_trade_name(chem_name):
    wiki_url = "https://en.wikipedia.org/wiki/"
    result = requests.get(wiki_url + chem_name.replace(" ", "+"))
    soup = bs(result.content, "html.parser")
    name = soup.find("h1", {"class": "firstHeading"}).text
    print(name)


def get_prices(chem_name):
    # Sigma
    sigma_url = "https://www.sigmaaldrich.com/catalog/search?term={}&interface=All".format(chem_name.replace(" ", "+"))
    print(sigma_url)
    # Right Price Chemicals
    right_price_url = "https://www.rightpricechemicals.com/catalogsearch/result/?q={}".format(
        chem_name.replace(" ", "+"))
    print(right_price_url)
    # Amazon
    amazon_url = "https://www.amazon.com/s?k={}&rh=p_72%3A2661619011".format(chem_name.replace(" ", "+"))
    print(amazon_url)
    # EBay
    ebay_url = "https://www.ebay.com/sch/i.html?_nkw={}".format(chem_name.replace(" ", "+"))
    print(ebay_url)


def get_safety_sheets(chem_name):
    # PubChem
    pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/#query={}".format(chem_name.replace(" ", "%20"))
    pubchem_url = get_rendered(pubchem_url, (By.LINK_TEXT, 'Summary')) \
        .find("div", {"id": "featured-results"}).find("a")["href"]
    print(pubchem_url)
    # Wikipedia
    wiki_url = "https://en.wikipedia.org/wiki/{}".format(chem_name.replace(" ", "%20"))
    print(wiki_url)
    # Sigma
    sigma_url = "https://www.sigmaaldrich.com/catalog/search?term={}&interface=All".format(chem_name.replace(" ", "+"))
    sigma_html = get_rendered(sigma_url, (By.CLASS_NAME, "msdsBulletPoint"))
    sigma_msds_id = sigma_html.find("a", {"class": "msdsBulletPoint"})["href"]\
        .split("(")[1][:-1].replace("'", "").split(",")[2][1:]
    sigma_url = "https://www.sigmaaldrich.com/MSDS/MSDS/DisplayMSDSPage.do?country=US&language=en&productNumber={}&brand=SIGALD".format(sigma_msds_id)
    print(sigma_url)


