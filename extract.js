const consideredNaringsvarden = ["Energi (kJ)", "Fett", "Summa mättade fettsyror", "Kolhydrater", "Socker totalt", "Fibrer", "Protein", "Salt"];

async function parseLivsmedelFromXML(filename) {
  const parser = new DOMParser();
  const data = await readFile(filename);
  const xmlDoc = parser.parseFromString(data, "text/xml");

  const result = [];
  for (const entry of xmlDoc.getElementsByTagName("Livsmedel")) {
    const livsmedelsnamn = entry.getElementsByTagName("Namn")[0].textContent;
    const naringsvarden = {};
    for (const naringsvarde of entry.getElementsByTagName("Naringsvarde")) {
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
    result.push({ livsmedelsnamn, naringsvarden });
  }
  result.sort(function(a, b) { return a.livsmedelsnamn.localeCompare(b.livsmedelsnamn, 'sv'); });

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

  for (const { livsmedelsnamn, naringsvarden } of livsmedel) {
    const livsmedelElement = createChildElement(xmlDoc, xmlDoc.documentElement, 'Livsmedel', '');
    createChildElement(xmlDoc, livsmedelElement, 'Namn', livsmedelsnamn);

    for (const naringsvarde in naringsvarden) {
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

async function populateFoods() {
  const livsmedel = await parseLivsmedelFromXML('compressed');
  const foodsSelect = document.querySelector("#foods");
  foodsSelect.innerHTML = "";
  for (const { livsmedelsnamn, naringsvarden } of livsmedel) {
    const option = document.createElement("option");
    option.textContent = livsmedelsnamn;
    const mappedArray = Object.entries(naringsvarden).map(([key, value]) => `${key}, ${value[0]}`);
    option.setAttribute('data-value', mappedArray.join(", "));
    foodsSelect.appendChild(option);
  }
}

function addFoods() {
  const selectedFoods = [];
  const foodsSelect = document.querySelector("#foods");
  const selectedOptions = foodsSelect.selectedOptions;
  for (const option of selectedOptions) {
    const foodName = option.textContent;
    selectedFoods.push([foodName, option.dataset.value]);
  }
  const selectedFoodsList = document.querySelector("#selected-foods");
  for (const food of selectedFoods) {
    const option = document.createElement("option");
    option.textContent = food[0];
    option.setAttribute('data-value', food[1]);
    selectedFoodsList.appendChild(option);
  }
}

function removeFoods() {
  const select = document.querySelector("#selected-foods");
  for (var i = select.options.length - 1; i >= 0; i--) {
    if (select.options[i].selected) {
      select.remove(i);
    }
  }
}

function searchFoods() {
  const input = document.getElementById('searchInput');
  const filter = input.value.toUpperCase();
  const foodsList = document.querySelector('#foods');
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
document.getElementById("add-food").addEventListener("click", addFoods);
document.getElementById("remove-food").addEventListener("click", removeFoods);
window.onload = populateFoods;
