# placeholder, still figuring this out

# parse url
library(XML)
library(RCurl)
library(rvest)
# parse url
url_parsed <- htmlParse(getURL("https://mycocosm.jgi.doe.gov/pages/fungi-1000-projects.jsf"), asText = TRUE)
tableNodes <- getNodeSet(url_parsed, "/html/body/div[3]/form/table[4]")
in_table <- readHTMLTable(tableNodes[[1]])
head(in_table)


innerTableNodes <- getNodeSet(url_parsed, "/html/body/div[3]/form/table[4]/tbody/tr[1]/td[6]/table")
inner_table <- readHTMLTable(innerTableNodes[[1]])


innerTableNodes[[1]] %>% html_nodes( "appropriate xpath or selector") %>% html_attr("href");

linkNodes <- getNodeSet(url_parsed, '//*/a')
pg %>% html_nodes("a") %>% html_attr("href")
linkNodes <- getNodeSet(url_parsed, "a")
in_link <- readHTMLTable(linkNodes[[1]])


webpage <- read_html("https://mycocosm.jgi.doe.gov/pages/fungi-1000-projects.jsf")
tbl_node <- webpage %>%
	html_nodes("table") %>%
	.[6]
link_nodes <- webpage %>%
	html_nodes("tr")


