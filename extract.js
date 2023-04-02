const consideredNaringsvarden = ["Energi (kJ)", "Fett", "Summa mättade fettsyror", "Kolhydrater", "Socker totalt", "Fibrer", "Protein", "Salt"];

async function parseLivsmedelFromXML(filename) {
  const parser = new DOMParser();
  const data = await readFile(filename);
  const xmlDoc = parser.parseFromString(data, "text/xml");

  const result = {};
  for (let entry of xmlDoc.getElementsByTagName("Livsmedel")) {
    const livsmedelsnamn = entry.getElementsByTagName("Namn")[0].textContent;
    const naringsvarden = {};
    for (let naringsvarde of entry.getElementsByTagName("Naringsvarde")) {
      const namn = naringsvarde.getElementsByTagName("Namn")[0].textContent;
      const varde = naringsvarde.getElementsByTagName("Varde")[0].textContent;
      const enhet = naringsvarde.getElementsByTagName("Enhet")[0].textContent;
      if (consideredNaringsvarden.includes(namn)) {
        if (enhet !== (namn === "Energi (kJ)" ? "kJ" : "g")) {
          throw new Error(`Unknown unit: ${enhet}`);
        }
        naringsvarden[namn] = [parseFloat(varde), enhet];
      }
    }
    result[livsmedelsnamn] = naringsvarden;
  }

  console.log(result);
  return result;
}

async function readFile(filename) {
  const response = await fetch(`./data/${filename}`);
  if (!response.ok) {
    throw new Error('Could not read file');
  }
  const data = await response.text();
  return data;
}

function createChildElement(xmlDoc, parentElement, name, value) {
  const element = xmlDoc.createElement(name);
  element.textContent = value;
  parentElement.appendChild(element);
  return element;
}

async function getLivsmedelAsString(filename) {
  const livsmedel = await parseLivsmedelFromXML(filename);

  // create a new XML document
  const xmlDoc = new DOMParser().parseFromString('<?xml version="1.0" encoding="UTF-8"?><root></root>', 'application/xml');

  for (let livsmedelsnamn in livsmedel) {
    const livsmedelElement = createChildElement(xmlDoc, xmlDoc.documentElement, 'Livsmedel', '');
    createChildElement(xmlDoc, livsmedelElement, 'Namn', livsmedelsnamn);

    const naringsvarden = livsmedel[livsmedelsnamn];
    for (let naringsvarde in naringsvarden) {
      const naringsvardeElement = createChildElement(xmlDoc, livsmedelElement, 'Naringsvarde', '');
      createChildElement(xmlDoc, naringsvardeElement, 'Namn', naringsvarde);
      createChildElement(xmlDoc, naringsvardeElement, 'Varde', naringsvarden[naringsvarde][0]);
      createChildElement(xmlDoc, naringsvardeElement, 'Enhet', naringsvarden[naringsvarde][1]);
    }
  }

  // convert the XML document to a string
  const xmlString = new XMLSerializer().serializeToString(xmlDoc);
  console.log(xmlString);
  return xmlString;
}

function downloadStringAsFile(stringData, fileName) {
  const blob = new Blob([stringData], {type: 'text/plain'});
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.download = fileName;
  link.href = url;
  link.click();
  URL.revokeObjectURL(url);
}

async function convertDatabase() {
  try {
    const xmlString = await getLivsmedelAsString('compressed');
    //downloadStringAsFile(xmlString, 'test.txt');
  }
  catch (error) {
    console.log(`Error ${error.lineNumber}: ${error.message}`);
  }
}

function searchFoods() {
  const input = document.getElementById('searchInput');
  const filter = input.value.toUpperCase();
  const foodsList = document.getElementById('foodsList');
  const foods = foodsList.getElementsByTagName('option');

  for (let i = 0; i < foods.length; i++) {
    const name = foods[i].textContent || foods[i].innerText;
    if (name.toUpperCase().indexOf(filter) > -1) {
      foods[i].style.display = '';
    } else {
      foods[i].style.display = 'none';
    }
  }
}

document.getElementById("extract").addEventListener("click", convertDatabase);
