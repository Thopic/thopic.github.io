/* tslint:disable */
/* eslint-disable */
/**
* @param {string} url
* @returns {Promise<any>}
*/
export function fetch_file(url: string): Promise<any>;
/**
*/
export class ResultJs {
  free(): void;
/**
*/
  readonly aa_cdr3: string;
/**
*/
  readonly aligned_j: string;
/**
*/
  readonly aligned_v: string;
/**
*/
  readonly cdr3_pos_string: string;
/**
*/
  readonly full_seq: string;
/**
*/
  readonly j_name: string;
/**
*/
  readonly likelihood: number;
/**
*/
  readonly n_cdr3: string;
/**
*/
  readonly pgen: number;
/**
*/
  readonly pos_start_cdr3: number;
/**
*/
  readonly reconstructed_seq: string;
/**
*/
  readonly seq: string;
/**
*/
  readonly v_name: string;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
  readonly memory: WebAssembly.Memory;
  readonly __wbg_resultjs_free: (a: number) => void;
  readonly resultjs_n_cdr3: (a: number, b: number) => void;
  readonly resultjs_aa_cdr3: (a: number, b: number) => void;
  readonly resultjs_seq: (a: number, b: number) => void;
  readonly resultjs_full_seq: (a: number, b: number) => void;
  readonly resultjs_reconstructed_seq: (a: number, b: number) => void;
  readonly resultjs_aligned_v: (a: number, b: number) => void;
  readonly resultjs_aligned_j: (a: number, b: number) => void;
  readonly resultjs_v_name: (a: number, b: number) => void;
  readonly resultjs_j_name: (a: number, b: number) => void;
  readonly resultjs_pgen: (a: number) => number;
  readonly resultjs_likelihood: (a: number) => number;
  readonly resultjs_pos_start_cdr3: (a: number) => number;
  readonly resultjs_cdr3_pos_string: (a: number, b: number) => void;
  readonly wrappedmodel_initialize: () => number;
  readonly wrappedmodel_load_model: (a: number, b: number, c: number, d: number) => number;
  readonly wrappedmodel_generate: (a: number, b: number) => void;
  readonly wrappedmodel_evaluate: (a: number, b: number, c: number) => void;
  readonly fetch_file: (a: number, b: number) => number;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
  readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
  readonly __wbindgen_export_2: WebAssembly.Table;
  readonly wasm_bindgen__convert__closures__invoke1_mut__h31b799c0ba683a3f: (a: number, b: number, c: number) => void;
  readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
  readonly __wbindgen_free: (a: number, b: number, c: number) => void;
  readonly __wbindgen_exn_store: (a: number) => void;
  readonly wasm_bindgen__convert__closures__invoke2_mut__h0f6d768d0ffb232d: (a: number, b: number, c: number, d: number) => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;
/**
* Instantiates the given `module`, which can either be bytes or
* a precompiled `WebAssembly.Module`.
*
* @param {SyncInitInput} module
*
* @returns {InitOutput}
*/
export function initSync(module: SyncInitInput): InitOutput;

/**
* If `module_or_path` is {RequestInfo} or {URL}, makes a request and
* for everything else, calls `WebAssembly.instantiate` directly.
*
* @param {InitInput | Promise<InitInput>} module_or_path
*
* @returns {Promise<InitOutput>}
*/
export default function __wbg_init (module_or_path?: InitInput | Promise<InitInput>): Promise<InitOutput>;
