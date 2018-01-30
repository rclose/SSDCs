#### DEFINE CONTINENTAL REGIONS ####
#Russia is currently included in Asia
#Cuba and Greenland are excluded from North America
regions <- list(
	c("United States","Canada","Mexico"),
	c("Argentina","Chile","Brazil","Bolivia","Colombia","Uruguay","Peru","Venezuela"),
	c("China","Mongolia","South Korea","Russian Federation","North Korea"),
	c("United Kingdom","France","Germany","Italy","Switzerland","Spain","Belgium","Germany","Romania","Sweden","Czech Republic","Denmark","Slovenia","Norway","Luxembourg","Netherlands","Ukraine","Hungary","Austria","Poland","Croatia","Portugal","Greece","Slovakia","Moldova, Republic of", "Serbia","Georgia","Ireland","Estonia"),
	c("Zambia","Namibia","Zimbabwe","Mali","Angola","Ethiopia","Cameroon","Malawi","Senegal","Tanzania","Eritrea","Sudan","Kenya","Libya","Niger","Tunisia","Algeria","Lesotho","Morocco","South Africa","Somalia","Djibouti","Gabon","Swaziland","Cote d'Ivoire","Mozambique","Congo, the Democratic Republic of the","Western Sahara","Nigeria","Ghana","Guinea","Madagascar")
)
bg_regions <- c("black","grey","brown","green","white")
names(regions) <- names(bg_regions) <- c("NAm","SA","AS","EU","AF")

library(countrycode)
regions <- lapply(regions, function(x) countrycode(x, origin = "country.name", destination = "iso2c"))
regions[[4]][1] <- "UK" #countrycode function uses GB instead of UK
