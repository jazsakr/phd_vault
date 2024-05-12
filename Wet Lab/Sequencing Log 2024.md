Create [7 color Gradient Palette here](https://coolors.co/gradient-palette/1f0021-751006?number=7).

![[sequencing_log_2024_heatmapcalender.png]]

```dataviewjs
dv.span("FSHD ;RiRunFill; ")

const searchString = "fshd"; // Search string to match

const calendarData = {
    year: 2024,
    colors: {
        color: ["#D5F1C8","#AAD496","#7FB764","#539A32","#287D00","#FF5800","#FF0000"], // w red for weekend
        //color: ["#D5F1C8","#B8DEA7","#9BCA85","#7FB764","#62A443","#459021","#287D00"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

// filter pages for "project" in yaml and matching string
for (let page of dv.pages('"Vault/daily note"').where(p => p.sequencing && p.day_of_week && p.sequencing.includes(searchString))) {
    // dv.span("page.file.name: "+page.file.name)
    calendarData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week
    });
}

// console.log(calendarData)

renderHeatmapCalendar(this.container, calendarData);
```

```dataviewjs
dv.span("FSHD ;RiRunFill; ")

const searchString = "fshd"; // Search string to match

const calendarData = {
    year: 2024,
    colors: {
        //color: ["#D5F1C8","#AAD496","#7FB764","#539A32","#287D00","#FF5800","#FF0000"], // w red for weekend
        color: ["#D5F1C8","#B8DEA7","#9BCA85","#7FB764","#62A443","#459021","#287D00"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

// filter pages for "project" in yaml and matching string
for (let page of dv.pages('"Vault/daily note"').where(p => p.sequencing && p.day_of_week && p.sequencing.includes(searchString))) {
    // dv.span("page.file.name: "+page.file.name)
    calendarData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week
    });
}

// console.log(calendarData)

renderHeatmapCalendar(this.container, calendarData);
```

```dataviewjs
dv.span("IGVF 🐁")

const searchString = "igvf"; // Search string to match

const calendarData = {
    year: 2024,
    colors: {
        color: ["#C7E8FF","#95CEF6","#64B4ED","#3299E3","#007FDA","#FF5800","#FF0000"], // w red for weekend
        //color: ["#C7E8FF","#A6D7F9","#85C5F3","#64B4ED","#42A2E6","#2191E0","#007FDA"],
    },
	intensityScaleStart: 1,
	intensityScaleEnd: 7,
    entries: []
}

// filter pages for "project" in yaml and matching string
for (let page of dv.pages('"Vault/daily note"').where(p => p.sequencing && p.day_of_week && p.sequencing.includes(searchString))) {
    // dv.span("page.file.name: "+page.file.name)
    calendarData.entries.push({
        date: page.file.name,
        intensity: page.day_of_week
    });
}

// console.log(calendarData)

renderHeatmapCalendar(this.container, calendarData);
```