const consideredNaringsvarden = ["Energi (kJ)", "Fett", "Summa mättade fettsyror", "Kolhydrater", "Socker totalt", "Fibrer", "Protein", "Salt"];

async function extract() {
  // Load the XML file
  const parser = new DOMParser();
  const data = await get_data();
  const xmlDoc = parser.parseFromString(data, "text/xml");

  // Loop through the entries and extract the values
  const result = [];
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
        naringsvarden[namn] = parseFloat(varde);
      }
    }
    result.push( { livsmedelsnamn, naringsvarden } );
  }

  console.log(result);
  return result;
}

async function get_data() {
  const response = await fetch('./data/20230401');
  if (!response.ok) {
    throw new Error('Could not read file');
  }
  const data = await response.text();
  return data;
}

function test() {
  try {
  }
  catch (error) {
    console.log(`Error ${error.lineNumber}: ${error.message}`);
  }
}

document.getElementById("extract").addEventListener("click", extract);
