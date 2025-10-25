```dataviewjs
dv.span("FSHD ;RiRunFill; ")

const searchString = "fshd"; // Search string to match

const trackerData = {
    separateMonths: true, 
    // heatmapTitle: "FSHD ;RiRunFill; ", 
    heatmapSubtitle: "Activity associated with the FSHD Project",
    colorScheme: {
        customColors: ["#D5F1C8","#AAD496","#7FB764","#539A32","#287D00","#FF5800","#FF0000"], // w red for weekend
        //color: ["#D5F1C8","#B8DEA7","#9BCA85","#7FB764","#62A443","#459021","#287D00"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

// filter pages for "project" in yaml and matching string
for (let page of dv.pages('"Vault/daily note"').where(p => p.project && p.day_of_week && p.project.includes(searchString))) {
    // dv.span("page.file.name: "+page.file.name)
    trackerData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week,
    }) 	
}

//console.log(trackerData)
	
renderHeatmapTracker(this.container, trackerData)
```

```dataviewjs
dv.span("IGVF 🐁")

const searchString = "igvf"; // Search string to match

const trackerData = {
    separateMonths: true, 
    // heatmapTitle: "IGVF 🐁", 
    heatmapSubtitle: "Activity associated with the IGVF Project",
    colorScheme: {
        customColors: ["#C7E8FF","#95CEF6","#64B4ED","#3299E3","#007FDA","#FF5800","#FF0000"], // w red for weekend
        //color: ["#C7E8FF","#A6D7F9","#85C5F3","#64B4ED","#42A2E6","#2191E0","#007FDA"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

// filter pages for "project" in yaml and matching string
for (let page of dv.pages('"Vault/daily note"').where(p => p.project && p.day_of_week && p.project.includes(searchString))) {
    // dv.span("page.file.name: "+page.file.name)
    trackerData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week,
    }) 	
}

//console.log(trackerData)
	
renderHeatmapTracker(this.container, trackerData)
```