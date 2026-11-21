/* Biowasm/Aioli tool wrappers. Aioli is loaded from CDN in index.html. */

let cli = null;
const requestedTools = new Set();

function parseSpec(spec) {
  const p = spec.split('/');
  return p.length === 3
    ? { tool: p[0], version: p[2], program: p[1], loading: 'lazy' }
    : { tool: p[0], version: p[1], loading: 'lazy' };
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
  cli = await new Aioli([...requestedTools].map(parseSpec), { printInterleaved: false });
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
    throw new Error(`aligner produced no output (stdout: ${(result.stdout || '').slice(-200)})`);
  }
  return out;
}

/** Run local BLAST against a user-provided FASTA database. */
export async function runBlast(program, queryFasta, dbFasta, maxHits) {
  await initTools(['blast/2.17.0']);
  await cli.mount([
    { name: 'query.fa', data: queryFasta },
    { name: 'db.fa', data: dbFasta },
  ]);
  // Build a local BLAST database from the mounted FASTA
  const dbType = program === 'blastn' ? 'nucl' : 'prot';
  await cli.exec(`makeblastdb -in db.fa -dbtype ${dbType} -out mydb`);
  const result = await cli.exec(
    `${program} -query query.fa -db mydb -max_target_seqs ${maxHits} -outfmt 6`
  );
  return result.stdout || result;
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
  return result.stdout || result;
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
