```dataviewjs
dv.span("Wet Lab 🥼")

const trackerData = {
	//year: 2025,
	separateMonths: true,
	//heatmapTitle: "Wet Lab 🥼",
	heatmapSubtitle: "Wet lab for all experiments",
	colorScheme: {
		customColors: ["#FF8800","#CFA22C","#9EBB58","#6ED584","#3DEEB0","#FF0000","#FF0000"], // with red for weekend
	},
	intensityScaleStart: 0,
	intensityScaleEnd: 6,   
	entries: []
}

for(let page of dv.pages('"Vault/Daily Note"')){

	// Get date and day of week from filename
	const dateStr = page.file.name;
	const date = new Date(dateStr);
	let dayOfWeek = date.getDay(); // (0 = Monday, ..., 5 = Saturday, 6 = Sunday)

	// Set icon per condition. Only one icon per day, in priority order.
	const hasSick = page.sick_day === true;
	const hasDayOff = page.day_off === true;
	const hasHoliday = page.holiday === true || typeof page.holiday === 'string';
	const hasConf = page.conference === true || typeof page.conference === 'string';
	
	let content = "";
	
	// Set content based on conditions
	if (hasHoliday) {content = "❎";
	} else if (hasDayOff) {content = "❌";
	} else if (hasSick) {content = "🤢";
	} else if (hasConf) {content = "🛩️";
	}
	
	// Create entry with proper structure
	const entry = {
		date: page.file.name,
		intensity: dayOfWeek,
		content: content
	};
	
	// Apply color scheme or dark gray default color
	if (page.wet_lab) {
		entry.customColor = undefined;
	} else {
		entry.customColor = "#333333";
	}
	
	trackerData.entries.push(entry);
}

renderHeatmapTracker(this.container, trackerData)
```

---

```dataviewjs
dv.span("Dry Lab 🖥")

const trackerData = {
	//year: 2025,
	separateMonths: true,
	//heatmapTitle: "Dry Lab 🖥",
	heatmapSubtitle: "Dry lab for all experiments",
	colorScheme: {
		customColors: ["#AFF2D8","#A6D7E2","#9CBCEC","#93A0F5","#8985FF", "#FF0000", "#FF0000"], // with red for weekend
	},
	intensityScaleStart: 0,
	intensityScaleEnd: 6,
	entries: []
}

for(let page of dv.pages('"Vault/Daily Note"')){

	// Get date and day of week from filename
	const dateStr = page.file.name;
	const date = new Date(dateStr);
	let dayOfWeek = date.getDay(); // (0 = Monday, ..., 5 = Saturday, 6 = Sunday)

	// Set icon per condition. Only one icon per day, in priority order.
	const hasSick = page.sick_day === true;
	const hasDayOff = page.day_off === true;
	const hasHoliday = page.holiday === true || typeof page.holiday === 'string';
	const hasConf = page.conference === true || typeof page.conference === 'string';
	
	let content = "";
	
	// Set content based on conditions
	if (hasHoliday) {content = "❎";
	} else if (hasDayOff) {content = "❌";
	} else if (hasSick) {content = "🤢";
	} else if (hasConf) {content = "🛩️";
	}
	
	// Create entry with proper structure
	const entry = {
		date: page.file.name,
		intensity: dayOfWeek,
		content: content
	};
	
	// Apply color scheme or dark gray default color
	if (page.dry_lab) {
		entry.customColor = undefined;
	} else {
		entry.customColor = "#333333";
	}
	
	trackerData.entries.push(entry);
}

renderHeatmapTracker(this.container, trackerData)
```

---

```dataviewjs
dv.span("Sequencing 🧬")

const trackerData = {
	//year: 2025,
	separateMonths: true,
	//heatmapTitle: "Sequencing 🧬",
	heatmapSubtitle: "Nanopore sequencing for all experiments",
	colorScheme: {
		customColors: ["#FFF34D","#CED664","#9CBA7B","#6B9D91","#3980A8", "#FF0000", "#FF0000"], // with red for weekend
	},
	intensityScaleStart: 0,
	intensityScaleEnd: 6,
	entries: []
}

for(let page of dv.pages('"Vault/Daily Note"')){

	// Get date and day of week from filename
	const dateStr = page.file.name;
	const date = new Date(dateStr);
	let dayOfWeek = date.getDay(); // (0 = Monday, ..., 5 = Saturday, 6 = Sunday)

	// Set icon per condition. Only one icon per day, in priority order.
	const hasSick = page.sick_day === true;
	const hasDayOff = page.day_off === true;
	const hasHoliday = page.holiday === true || typeof page.holiday === 'string';
	const hasConf = page.conference === true || typeof page.conference === 'string';
	
	let content = "";
	
	if (hasHoliday) {content = "❎";
	} else if (hasDayOff) {content = "❌";
	} else if (hasSick) {content = "🤢";
	} else if (hasConf) {content = "🛩️";
	}
	
	// Create entry with proper structure
	const entry = {
		date: page.file.name,
		intensity: dayOfWeek,
		content: content
	};
	
	// Apply color scheme or dark gray default color
	if (page.sequencing) {
		entry.customColor = undefined;
	} else {
		entry.customColor = "#333333";
	}
	
	trackerData.entries.push(entry);
}

renderHeatmapTracker(this.container, trackerData)
```