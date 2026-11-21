/* Biowasm/Aioli tool wrappers. Aioli is loaded from CDN in index.html. */

let cli = null;
const requestedTools = new Set();

function parseSpec(spec) {
  const p = spec.split('/');
  return p.length === 3
    ? { tool: p[0], version: p[2], program: p[1] }
    : { tool: p[0], version: p[1] };
}

/** Normalize Aioli exec() output (plain string, or {stdout, stderr} object). */
function execText(result) {
  if (typeof result === 'string') return result;
  return result?.stdout ?? '';
}

export async function initTools(tools) {
  if (typeof Aioli === 'undefined') {
    throw new Error('Aioli failed to load from CDN — check your network connection');
  }
  const before = requestedTools.size;
  tools.forEach(t => requestedTools.add(t));
  if (cli && requestedTools.size === before) return cli;
  // Aioli's filesystem is shared per instance; re-init with the full tool set
  // whenever a new tool is added so earlier tools remain available.
  // Use Aioli's default output mode: exec() returns a plain string. Some
  // programs (notably BLAST) do not populate stdout when printInterleaved=false.
  cli = await new Aioli([...requestedTools].map(parseSpec));
  return cli;
}

/** Mount a FASTA string as a virtual file, run an aligner, return aligned FASTA. */
export async function runAligner(tool, fasta) {
  const specs = {
    muscle: 'muscle/5.1.0',
    kalign: 'kalign/3.3.1',
  };
  await initTools([specs[tool]]);
  const inputName = 'input.fa';
  await cli.mount([{ name: inputName, data: fasta }]);

  let cmd;
  switch (tool) {
    case 'muscle':
      cmd = `muscle -align ${inputName} -output aligned.fa`;
      break;
    case 'kalign':
      cmd = `kalign -i ${inputName} -o aligned.fa`;
      break;
    default:
      throw new Error(`Unknown aligner: ${tool}`);
  }

  const result = await cli.exec(cmd);
  // Biowasm tools run with cwd=/shared/data; output files are relative to it.
  const out = await cli.cat('aligned.fa');
  if (typeof out !== 'string' || !out.trim()) {
    throw new Error(`aligner produced no output (stdout: ${execText(result).slice(-200)})`);
  }
  return out;
}

/** Run local BLAST against a user-provided FASTA database. */
export async function runBlast(program, queryFasta, dbFasta, maxHits) {
  // Biowasm exposes each BLAST binary as a separate program (there is no
  // monolithic `blast.js`), so request the chosen program plus makeblastdb.
  await initTools([`blast/${program}/2.17.0`, 'blast/makeblastdb/2.17.0']);
  await cli.mount([
    { name: 'query.fa', data: queryFasta },
    { name: 'db.fa', data: dbFasta },
  ]);
  // Build a local BLAST database from the mounted FASTA
  const isNucleotide = program === 'blastn';
  const dbType = isNucleotide ? 'nucl' : 'prot';
  // v5 databases use LMDB, which is unavailable in the WASM filesystem;
  // v4 uses traditional .pin/.psq index files and works in-browser.
  await cli.exec(`makeblastdb -in db.fa -dbtype ${dbType} -out mydb -blastdb_version 4`);
  // Use a small word size for short queries; disable low-complexity masking so
  // tandem-repeat queries aren't filtered out entirely.
  const filter = isNucleotide ? '-task blastn -dust no' : '-seg no';
  const result = await cli.exec(
    `${program} -query query.fa -db mydb ${filter} -max_target_seqs ${maxHits} -outfmt 6`
  );
  const text = execText(result);
  if (!text.trim()) throw new Error(`${program} produced no hits — try a different query/database`);
  return text;
}

/** Run local phmmer against a user-provided FASTA database. */
export async function runPhmmer(queryFasta, dbFasta, maxHits) {
  await initTools(['hmmer/3.4']);
  await cli.mount([
    { name: 'query.fa', data: queryFasta },
    { name: 'db.fa', data: dbFasta },
  ]);
  const result = await cli.exec(
    `phmmer --maxhits ${maxHits} --domT 1e-3 query.fa db.fa`
  );
  const text = execText(result);
  if (!text.trim()) throw new Error('phmmer produced no hits — try a different query/database');
  return text;
}

/** Parse BLAST tabular output -> array of subject IDs. */
export function parseBlastHits(tabular) {
  const ids = [];
  for (const line of tabular.split('\n')) {
    if (!line.trim() || line.startsWith('#')) continue;
    const fields = line.split('\t');
    if (fields.length > 1) ids.push(fields[1]);
  }
  return [...new Set(ids)];
}

/** Parse phmmer output -> array of subject names (lines starting with '>>'). */
export function parsePhmmerHits(output) {
  const ids = [];
  for (const line of output.split('\n')) {
    if (line.startsWith('>>')) ids.push(line.slice(2).trim().split(/\s+/)[0]);
  }
  return [...new Set(ids)];
}
