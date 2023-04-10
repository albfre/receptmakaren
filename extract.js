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
    const xmlString = await getLivsmedelAsString('20230401');
    downloadStringAsFile(xmlString, 'compressed');
  }
  catch (error) {
    console.log(`Error ${error.lineNumber}: ${error.message}`);
  }
}

document.getElementById("extract").addEventListener("click", convertDatabase);
