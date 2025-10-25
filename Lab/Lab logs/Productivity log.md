```dataviewjs
dv.span("Sequencing 🧬")

const trackerData = {
    separateMonths: true, 
    // heatmapTitle: "Sequencing 🧬", 
    heatmapSubtitle: "Days Nanopore sequencing for all experiments",
    colorScheme: {
        customColors: ["#FFF34D","#CED664","#9CBA7B","#6B9D91","#3980A8", "#FF0000", "#FF0000"], // w red for weekend
        // ["#FFF34D","#DEE05C","#BDCD6B","#9CBA7B","#7BA68A", "#5A9399", "#3980A8"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

for(let page of dv.pages('"Vault/daily note"').where(p=>p.sequencing)){
	
    trackerData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week,
    }) 	
}

//console.log(trackerData)
	
renderHeatmapTracker(this.container, trackerData)

```




```dataviewjs
dv.span("Wet Lab 🥼")

const trackerData = {
    separateMonths: true, 
    // heatmapTitle: "Wet Lab 🥼", 
    heatmapSubtitle: "Wet lab for all experiments",
    colorScheme: {
        customColors: ["#FF8800","#CFA22C","#9EBB58","#6ED584","#3DEEB0", "#FF0000", "#FF0000"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

for(let page of dv.pages('"Vault/daily note"').where(p=>p.wet_lab && p.day_of_week)){
	 
    trackerData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week,
    }) 	
}

//console.log(trackerData)
	
renderHeatmapTracker(this.container, trackerData)

```



```dataviewjs
dv.span("Dry Lab 🖥")

const trackerData = {
    separateMonths: true, 
    // heatmapTitle: "Dry Lab 🖥", 
    heatmapSubtitle: "Dry lab for all experiments",
    colorScheme: {
        customColors: ["#AFF2D8","#A6D7E2","#9CBCEC","#93A0F5","#8985FF", "#FF0000", "#FF0000"]
        // greyscale: ["#DDDDDD","#BEBEBE","#A0A0A0","#818181","#626262", "#444444", "#252525"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

for(let page of dv.pages('"Vault/daily note"').where(p=>p.dry_lab && p.day_of_week)){
	 
    trackerData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week,
    }) 	
}

//console.log(trackerData)
	
renderHeatmapTracker(this.container, trackerData)

```