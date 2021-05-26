<!DOCTYPE html>
<html>
<head><meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0">

<title>Nellie_D17</title><script src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.1.10/require.min.js"></script>




<style type="text/css">
    pre { line-height: 125%; }
td.linenos .normal { color: inherit; background-color: transparent; padding-left: 5px; padding-right: 5px; }
span.linenos { color: inherit; background-color: transparent; padding-left: 5px; padding-right: 5px; }
td.linenos .special { color: #000000; background-color: #ffffc0; padding-left: 5px; padding-right: 5px; }
span.linenos.special { color: #000000; background-color: #ffffc0; padding-left: 5px; padding-right: 5px; }
.highlight .hll { background-color: var(--jp-cell-editor-active-background) }
.highlight { background: var(--jp-cell-editor-background); color: var(--jp-mirror-editor-variable-color) }
.highlight .c { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment */
.highlight .err { color: var(--jp-mirror-editor-error-color) } /* Error */
.highlight .k { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword */
.highlight .o { color: var(--jp-mirror-editor-operator-color); font-weight: bold } /* Operator */
.highlight .p { color: var(--jp-mirror-editor-punctuation-color) } /* Punctuation */
.highlight .ch { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment.Hashbang */
.highlight .cm { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment.Multiline */
.highlight .cp { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment.Preproc */
.highlight .cpf { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment.PreprocFile */
.highlight .c1 { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment.Single */
.highlight .cs { color: var(--jp-mirror-editor-comment-color); font-style: italic } /* Comment.Special */
.highlight .kc { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword.Constant */
.highlight .kd { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword.Declaration */
.highlight .kn { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword.Namespace */
.highlight .kp { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword.Pseudo */
.highlight .kr { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword.Reserved */
.highlight .kt { color: var(--jp-mirror-editor-keyword-color); font-weight: bold } /* Keyword.Type */
.highlight .m { color: var(--jp-mirror-editor-number-color) } /* Literal.Number */
.highlight .s { color: var(--jp-mirror-editor-string-color) } /* Literal.String */
.highlight .ow { color: var(--jp-mirror-editor-operator-color); font-weight: bold } /* Operator.Word */
.highlight .w { color: var(--jp-mirror-editor-variable-color) } /* Text.Whitespace */
.highlight .mb { color: var(--jp-mirror-editor-number-color) } /* Literal.Number.Bin */
.highlight .mf { color: var(--jp-mirror-editor-number-color) } /* Literal.Number.Float */
.highlight .mh { color: var(--jp-mirror-editor-number-color) } /* Literal.Number.Hex */
.highlight .mi { color: var(--jp-mirror-editor-number-color) } /* Literal.Number.Integer */
.highlight .mo { color: var(--jp-mirror-editor-number-color) } /* Literal.Number.Oct */
.highlight .sa { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Affix */
.highlight .sb { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Backtick */
.highlight .sc { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Char */
.highlight .dl { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Delimiter */
.highlight .sd { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Doc */
.highlight .s2 { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Double */
.highlight .se { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Escape */
.highlight .sh { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Heredoc */
.highlight .si { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Interpol */
.highlight .sx { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Other */
.highlight .sr { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Regex */
.highlight .s1 { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Single */
.highlight .ss { color: var(--jp-mirror-editor-string-color) } /* Literal.String.Symbol */
.highlight .il { color: var(--jp-mirror-editor-number-color) } /* Literal.Number.Integer.Long */
  </style>



<style type="text/css">
/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*
 * Mozilla scrollbar styling
 */

/* use standard opaque scrollbars for most nodes */
[data-jp-theme-scrollbars='true'] {
  scrollbar-color: rgb(var(--jp-scrollbar-thumb-color))
    var(--jp-scrollbar-background-color);
}

/* for code nodes, use a transparent style of scrollbar. These selectors
 * will match lower in the tree, and so will override the above */
[data-jp-theme-scrollbars='true'] .CodeMirror-hscrollbar,
[data-jp-theme-scrollbars='true'] .CodeMirror-vscrollbar {
  scrollbar-color: rgba(var(--jp-scrollbar-thumb-color), 0.5) transparent;
}

/*
 * Webkit scrollbar styling
 */

/* use standard opaque scrollbars for most nodes */

[data-jp-theme-scrollbars='true'] ::-webkit-scrollbar,
[data-jp-theme-scrollbars='true'] ::-webkit-scrollbar-corner {
  background: var(--jp-scrollbar-background-color);
}

[data-jp-theme-scrollbars='true'] ::-webkit-scrollbar-thumb {
  background: rgb(var(--jp-scrollbar-thumb-color));
  border: var(--jp-scrollbar-thumb-margin) solid transparent;
  background-clip: content-box;
  border-radius: var(--jp-scrollbar-thumb-radius);
}

[data-jp-theme-scrollbars='true'] ::-webkit-scrollbar-track:horizontal {
  border-left: var(--jp-scrollbar-endpad) solid
    var(--jp-scrollbar-background-color);
  border-right: var(--jp-scrollbar-endpad) solid
    var(--jp-scrollbar-background-color);
}

[data-jp-theme-scrollbars='true'] ::-webkit-scrollbar-track:vertical {
  border-top: var(--jp-scrollbar-endpad) solid
    var(--jp-scrollbar-background-color);
  border-bottom: var(--jp-scrollbar-endpad) solid
    var(--jp-scrollbar-background-color);
}

/* for code nodes, use a transparent style of scrollbar */

[data-jp-theme-scrollbars='true'] .CodeMirror-hscrollbar::-webkit-scrollbar,
[data-jp-theme-scrollbars='true'] .CodeMirror-vscrollbar::-webkit-scrollbar,
[data-jp-theme-scrollbars='true']
  .CodeMirror-hscrollbar::-webkit-scrollbar-corner,
[data-jp-theme-scrollbars='true']
  .CodeMirror-vscrollbar::-webkit-scrollbar-corner {
  background-color: transparent;
}

[data-jp-theme-scrollbars='true']
  .CodeMirror-hscrollbar::-webkit-scrollbar-thumb,
[data-jp-theme-scrollbars='true']
  .CodeMirror-vscrollbar::-webkit-scrollbar-thumb {
  background: rgba(var(--jp-scrollbar-thumb-color), 0.5);
  border: var(--jp-scrollbar-thumb-margin) solid transparent;
  background-clip: content-box;
  border-radius: var(--jp-scrollbar-thumb-radius);
}

[data-jp-theme-scrollbars='true']
  .CodeMirror-hscrollbar::-webkit-scrollbar-track:horizontal {
  border-left: var(--jp-scrollbar-endpad) solid transparent;
  border-right: var(--jp-scrollbar-endpad) solid transparent;
}

[data-jp-theme-scrollbars='true']
  .CodeMirror-vscrollbar::-webkit-scrollbar-track:vertical {
  border-top: var(--jp-scrollbar-endpad) solid transparent;
  border-bottom: var(--jp-scrollbar-endpad) solid transparent;
}

/*
 * Phosphor
 */

.lm-ScrollBar[data-orientation='horizontal'] {
  min-height: 16px;
  max-height: 16px;
  min-width: 45px;
  border-top: 1px solid #a0a0a0;
}

.lm-ScrollBar[data-orientation='vertical'] {
  min-width: 16px;
  max-width: 16px;
  min-height: 45px;
  border-left: 1px solid #a0a0a0;
}

.lm-ScrollBar-button {
  background-color: #f0f0f0;
  background-position: center center;
  min-height: 15px;
  max-height: 15px;
  min-width: 15px;
  max-width: 15px;
}

.lm-ScrollBar-button:hover {
  background-color: #dadada;
}

.lm-ScrollBar-button.lm-mod-active {
  background-color: #cdcdcd;
}

.lm-ScrollBar-track {
  background: #f0f0f0;
}

.lm-ScrollBar-thumb {
  background: #cdcdcd;
}

.lm-ScrollBar-thumb:hover {
  background: #bababa;
}

.lm-ScrollBar-thumb.lm-mod-active {
  background: #a0a0a0;
}

.lm-ScrollBar[data-orientation='horizontal'] .lm-ScrollBar-thumb {
  height: 100%;
  min-width: 15px;
  border-left: 1px solid #a0a0a0;
  border-right: 1px solid #a0a0a0;
}

.lm-ScrollBar[data-orientation='vertical'] .lm-ScrollBar-thumb {
  width: 100%;
  min-height: 15px;
  border-top: 1px solid #a0a0a0;
  border-bottom: 1px solid #a0a0a0;
}

.lm-ScrollBar[data-orientation='horizontal']
  .lm-ScrollBar-button[data-action='decrement'] {
  background-image: var(--jp-icon-caret-left);
  background-size: 17px;
}

.lm-ScrollBar[data-orientation='horizontal']
  .lm-ScrollBar-button[data-action='increment'] {
  background-image: var(--jp-icon-caret-right);
  background-size: 17px;
}

.lm-ScrollBar[data-orientation='vertical']
  .lm-ScrollBar-button[data-action='decrement'] {
  background-image: var(--jp-icon-caret-up);
  background-size: 17px;
}

.lm-ScrollBar[data-orientation='vertical']
  .lm-ScrollBar-button[data-action='increment'] {
  background-image: var(--jp-icon-caret-down);
  background-size: 17px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-Widget, /* </DEPRECATED> */
.lm-Widget {
  box-sizing: border-box;
  position: relative;
  overflow: hidden;
  cursor: default;
}


/* <DEPRECATED> */ .p-Widget.p-mod-hidden, /* </DEPRECATED> */
.lm-Widget.lm-mod-hidden {
  display: none !important;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-CommandPalette, /* </DEPRECATED> */
.lm-CommandPalette {
  display: flex;
  flex-direction: column;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}


/* <DEPRECATED> */ .p-CommandPalette-search, /* </DEPRECATED> */
.lm-CommandPalette-search {
  flex: 0 0 auto;
}


/* <DEPRECATED> */ .p-CommandPalette-content, /* </DEPRECATED> */
.lm-CommandPalette-content {
  flex: 1 1 auto;
  margin: 0;
  padding: 0;
  min-height: 0;
  overflow: auto;
  list-style-type: none;
}


/* <DEPRECATED> */ .p-CommandPalette-header, /* </DEPRECATED> */
.lm-CommandPalette-header {
  overflow: hidden;
  white-space: nowrap;
  text-overflow: ellipsis;
}


/* <DEPRECATED> */ .p-CommandPalette-item, /* </DEPRECATED> */
.lm-CommandPalette-item {
  display: flex;
  flex-direction: row;
}


/* <DEPRECATED> */ .p-CommandPalette-itemIcon, /* </DEPRECATED> */
.lm-CommandPalette-itemIcon {
  flex: 0 0 auto;
}


/* <DEPRECATED> */ .p-CommandPalette-itemContent, /* </DEPRECATED> */
.lm-CommandPalette-itemContent {
  flex: 1 1 auto;
  overflow: hidden;
}


/* <DEPRECATED> */ .p-CommandPalette-itemShortcut, /* </DEPRECATED> */
.lm-CommandPalette-itemShortcut {
  flex: 0 0 auto;
}


/* <DEPRECATED> */ .p-CommandPalette-itemLabel, /* </DEPRECATED> */
.lm-CommandPalette-itemLabel {
  overflow: hidden;
  white-space: nowrap;
  text-overflow: ellipsis;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-DockPanel, /* </DEPRECATED> */
.lm-DockPanel {
  z-index: 0;
}


/* <DEPRECATED> */ .p-DockPanel-widget, /* </DEPRECATED> */
.lm-DockPanel-widget {
  z-index: 0;
}


/* <DEPRECATED> */ .p-DockPanel-tabBar, /* </DEPRECATED> */
.lm-DockPanel-tabBar {
  z-index: 1;
}


/* <DEPRECATED> */ .p-DockPanel-handle, /* </DEPRECATED> */
.lm-DockPanel-handle {
  z-index: 2;
}


/* <DEPRECATED> */ .p-DockPanel-handle.p-mod-hidden, /* </DEPRECATED> */
.lm-DockPanel-handle.lm-mod-hidden {
  display: none !important;
}


/* <DEPRECATED> */ .p-DockPanel-handle:after, /* </DEPRECATED> */
.lm-DockPanel-handle:after {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  content: '';
}


/* <DEPRECATED> */
.p-DockPanel-handle[data-orientation='horizontal'],
/* </DEPRECATED> */
.lm-DockPanel-handle[data-orientation='horizontal'] {
  cursor: ew-resize;
}


/* <DEPRECATED> */
.p-DockPanel-handle[data-orientation='vertical'],
/* </DEPRECATED> */
.lm-DockPanel-handle[data-orientation='vertical'] {
  cursor: ns-resize;
}


/* <DEPRECATED> */
.p-DockPanel-handle[data-orientation='horizontal']:after,
/* </DEPRECATED> */
.lm-DockPanel-handle[data-orientation='horizontal']:after {
  left: 50%;
  min-width: 8px;
  transform: translateX(-50%);
}


/* <DEPRECATED> */
.p-DockPanel-handle[data-orientation='vertical']:after,
/* </DEPRECATED> */
.lm-DockPanel-handle[data-orientation='vertical']:after {
  top: 50%;
  min-height: 8px;
  transform: translateY(-50%);
}


/* <DEPRECATED> */ .p-DockPanel-overlay, /* </DEPRECATED> */
.lm-DockPanel-overlay {
  z-index: 3;
  box-sizing: border-box;
  pointer-events: none;
}


/* <DEPRECATED> */ .p-DockPanel-overlay.p-mod-hidden, /* </DEPRECATED> */
.lm-DockPanel-overlay.lm-mod-hidden {
  display: none !important;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-Menu, /* </DEPRECATED> */
.lm-Menu {
  z-index: 10000;
  position: absolute;
  white-space: nowrap;
  overflow-x: hidden;
  overflow-y: auto;
  outline: none;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}


/* <DEPRECATED> */ .p-Menu-content, /* </DEPRECATED> */
.lm-Menu-content {
  margin: 0;
  padding: 0;
  display: table;
  list-style-type: none;
}


/* <DEPRECATED> */ .p-Menu-item, /* </DEPRECATED> */
.lm-Menu-item {
  display: table-row;
}


/* <DEPRECATED> */
.p-Menu-item.p-mod-hidden,
.p-Menu-item.p-mod-collapsed,
/* </DEPRECATED> */
.lm-Menu-item.lm-mod-hidden,
.lm-Menu-item.lm-mod-collapsed {
  display: none !important;
}


/* <DEPRECATED> */
.p-Menu-itemIcon,
.p-Menu-itemSubmenuIcon,
/* </DEPRECATED> */
.lm-Menu-itemIcon,
.lm-Menu-itemSubmenuIcon {
  display: table-cell;
  text-align: center;
}


/* <DEPRECATED> */ .p-Menu-itemLabel, /* </DEPRECATED> */
.lm-Menu-itemLabel {
  display: table-cell;
  text-align: left;
}


/* <DEPRECATED> */ .p-Menu-itemShortcut, /* </DEPRECATED> */
.lm-Menu-itemShortcut {
  display: table-cell;
  text-align: right;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-MenuBar, /* </DEPRECATED> */
.lm-MenuBar {
  outline: none;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}


/* <DEPRECATED> */ .p-MenuBar-content, /* </DEPRECATED> */
.lm-MenuBar-content {
  margin: 0;
  padding: 0;
  display: flex;
  flex-direction: row;
  list-style-type: none;
}


/* <DEPRECATED> */ .p--MenuBar-item, /* </DEPRECATED> */
.lm-MenuBar-item {
  box-sizing: border-box;
}


/* <DEPRECATED> */
.p-MenuBar-itemIcon,
.p-MenuBar-itemLabel,
/* </DEPRECATED> */
.lm-MenuBar-itemIcon,
.lm-MenuBar-itemLabel {
  display: inline-block;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-ScrollBar, /* </DEPRECATED> */
.lm-ScrollBar {
  display: flex;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}


/* <DEPRECATED> */
.p-ScrollBar[data-orientation='horizontal'],
/* </DEPRECATED> */
.lm-ScrollBar[data-orientation='horizontal'] {
  flex-direction: row;
}


/* <DEPRECATED> */
.p-ScrollBar[data-orientation='vertical'],
/* </DEPRECATED> */
.lm-ScrollBar[data-orientation='vertical'] {
  flex-direction: column;
}


/* <DEPRECATED> */ .p-ScrollBar-button, /* </DEPRECATED> */
.lm-ScrollBar-button {
  box-sizing: border-box;
  flex: 0 0 auto;
}


/* <DEPRECATED> */ .p-ScrollBar-track, /* </DEPRECATED> */
.lm-ScrollBar-track {
  box-sizing: border-box;
  position: relative;
  overflow: hidden;
  flex: 1 1 auto;
}


/* <DEPRECATED> */ .p-ScrollBar-thumb, /* </DEPRECATED> */
.lm-ScrollBar-thumb {
  box-sizing: border-box;
  position: absolute;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-SplitPanel-child, /* </DEPRECATED> */
.lm-SplitPanel-child {
  z-index: 0;
}


/* <DEPRECATED> */ .p-SplitPanel-handle, /* </DEPRECATED> */
.lm-SplitPanel-handle {
  z-index: 1;
}


/* <DEPRECATED> */ .p-SplitPanel-handle.p-mod-hidden, /* </DEPRECATED> */
.lm-SplitPanel-handle.lm-mod-hidden {
  display: none !important;
}


/* <DEPRECATED> */ .p-SplitPanel-handle:after, /* </DEPRECATED> */
.lm-SplitPanel-handle:after {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  content: '';
}


/* <DEPRECATED> */
.p-SplitPanel[data-orientation='horizontal'] > .p-SplitPanel-handle,
/* </DEPRECATED> */
.lm-SplitPanel[data-orientation='horizontal'] > .lm-SplitPanel-handle {
  cursor: ew-resize;
}


/* <DEPRECATED> */
.p-SplitPanel[data-orientation='vertical'] > .p-SplitPanel-handle,
/* </DEPRECATED> */
.lm-SplitPanel[data-orientation='vertical'] > .lm-SplitPanel-handle {
  cursor: ns-resize;
}


/* <DEPRECATED> */
.p-SplitPanel[data-orientation='horizontal'] > .p-SplitPanel-handle:after,
/* </DEPRECATED> */
.lm-SplitPanel[data-orientation='horizontal'] > .lm-SplitPanel-handle:after {
  left: 50%;
  min-width: 8px;
  transform: translateX(-50%);
}


/* <DEPRECATED> */
.p-SplitPanel[data-orientation='vertical'] > .p-SplitPanel-handle:after,
/* </DEPRECATED> */
.lm-SplitPanel[data-orientation='vertical'] > .lm-SplitPanel-handle:after {
  top: 50%;
  min-height: 8px;
  transform: translateY(-50%);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-TabBar, /* </DEPRECATED> */
.lm-TabBar {
  display: flex;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}


/* <DEPRECATED> */ .p-TabBar[data-orientation='horizontal'], /* </DEPRECATED> */
.lm-TabBar[data-orientation='horizontal'] {
  flex-direction: row;
}


/* <DEPRECATED> */ .p-TabBar[data-orientation='vertical'], /* </DEPRECATED> */
.lm-TabBar[data-orientation='vertical'] {
  flex-direction: column;
}


/* <DEPRECATED> */ .p-TabBar-content, /* </DEPRECATED> */
.lm-TabBar-content {
  margin: 0;
  padding: 0;
  display: flex;
  flex: 1 1 auto;
  list-style-type: none;
}


/* <DEPRECATED> */
.p-TabBar[data-orientation='horizontal'] > .p-TabBar-content,
/* </DEPRECATED> */
.lm-TabBar[data-orientation='horizontal'] > .lm-TabBar-content {
  flex-direction: row;
}


/* <DEPRECATED> */
.p-TabBar[data-orientation='vertical'] > .p-TabBar-content,
/* </DEPRECATED> */
.lm-TabBar[data-orientation='vertical'] > .lm-TabBar-content {
  flex-direction: column;
}


/* <DEPRECATED> */ .p-TabBar-tab, /* </DEPRECATED> */
.lm-TabBar-tab {
  display: flex;
  flex-direction: row;
  box-sizing: border-box;
  overflow: hidden;
}


/* <DEPRECATED> */
.p-TabBar-tabIcon,
.p-TabBar-tabCloseIcon,
/* </DEPRECATED> */
.lm-TabBar-tabIcon,
.lm-TabBar-tabCloseIcon {
  flex: 0 0 auto;
}


/* <DEPRECATED> */ .p-TabBar-tabLabel, /* </DEPRECATED> */
.lm-TabBar-tabLabel {
  flex: 1 1 auto;
  overflow: hidden;
  white-space: nowrap;
}


/* <DEPRECATED> */ .p-TabBar-tab.p-mod-hidden, /* </DEPRECATED> */
.lm-TabBar-tab.lm-mod-hidden {
  display: none !important;
}


/* <DEPRECATED> */ .p-TabBar.p-mod-dragging .p-TabBar-tab, /* </DEPRECATED> */
.lm-TabBar.lm-mod-dragging .lm-TabBar-tab {
  position: relative;
}


/* <DEPRECATED> */
.p-TabBar.p-mod-dragging[data-orientation='horizontal'] .p-TabBar-tab,
/* </DEPRECATED> */
.lm-TabBar.lm-mod-dragging[data-orientation='horizontal'] .lm-TabBar-tab {
  left: 0;
  transition: left 150ms ease;
}


/* <DEPRECATED> */
.p-TabBar.p-mod-dragging[data-orientation='vertical'] .p-TabBar-tab,
/* </DEPRECATED> */
.lm-TabBar.lm-mod-dragging[data-orientation='vertical'] .lm-TabBar-tab {
  top: 0;
  transition: top 150ms ease;
}


/* <DEPRECATED> */
.p-TabBar.p-mod-dragging .p-TabBar-tab.p-mod-dragging
/* </DEPRECATED> */
.lm-TabBar.lm-mod-dragging .lm-TabBar-tab.lm-mod-dragging {
  transition: none;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ .p-TabPanel-tabBar, /* </DEPRECATED> */
.lm-TabPanel-tabBar {
  z-index: 1;
}


/* <DEPRECATED> */ .p-TabPanel-stackedPanel, /* </DEPRECATED> */
.lm-TabPanel-stackedPanel {
  z-index: 0;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/

@charset "UTF-8";
/*!

Copyright 2015-present Palantir Technologies, Inc. All rights reserved.
Licensed under the Apache License, Version 2.0.

*/
html{
  -webkit-box-sizing:border-box;
          box-sizing:border-box; }

*,
*::before,
*::after{
  -webkit-box-sizing:inherit;
          box-sizing:inherit; }

body{
  text-transform:none;
  line-height:1.28581;
  letter-spacing:0;
  font-size:14px;
  font-weight:400;
  color:#182026;
  font-family:-apple-system, "BlinkMacSystemFont", "Segoe UI", "Roboto", "Oxygen", "Ubuntu", "Cantarell", "Open Sans", "Helvetica Neue", "Icons16", sans-serif; }

p{
  margin-top:0;
  margin-bottom:10px; }

small{
  font-size:12px; }

strong{
  font-weight:600; }

::-moz-selection{
  background:rgba(125, 188, 255, 0.6); }

::selection{
  background:rgba(125, 188, 255, 0.6); }
.bp3-heading{
  color:#182026;
  font-weight:600;
  margin:0 0 10px;
  padding:0; }
  .bp3-dark .bp3-heading{
    color:#f5f8fa; }

h1.bp3-heading, .bp3-running-text h1{
  line-height:40px;
  font-size:36px; }

h2.bp3-heading, .bp3-running-text h2{
  line-height:32px;
  font-size:28px; }

h3.bp3-heading, .bp3-running-text h3{
  line-height:25px;
  font-size:22px; }

h4.bp3-heading, .bp3-running-text h4{
  line-height:21px;
  font-size:18px; }

h5.bp3-heading, .bp3-running-text h5{
  line-height:19px;
  font-size:16px; }

h6.bp3-heading, .bp3-running-text h6{
  line-height:16px;
  font-size:14px; }
.bp3-ui-text{
  text-transform:none;
  line-height:1.28581;
  letter-spacing:0;
  font-size:14px;
  font-weight:400; }

.bp3-monospace-text{
  text-transform:none;
  font-family:monospace; }

.bp3-text-muted{
  color:#5c7080; }
  .bp3-dark .bp3-text-muted{
    color:#a7b6c2; }

.bp3-text-disabled{
  color:rgba(92, 112, 128, 0.6); }
  .bp3-dark .bp3-text-disabled{
    color:rgba(167, 182, 194, 0.6); }

.bp3-text-overflow-ellipsis{
  overflow:hidden;
  text-overflow:ellipsis;
  white-space:nowrap;
  word-wrap:normal; }
.bp3-running-text{
  line-height:1.5;
  font-size:14px; }
  .bp3-running-text h1{
    color:#182026;
    font-weight:600;
    margin-top:40px;
    margin-bottom:20px; }
    .bp3-dark .bp3-running-text h1{
      color:#f5f8fa; }
  .bp3-running-text h2{
    color:#182026;
    font-weight:600;
    margin-top:40px;
    margin-bottom:20px; }
    .bp3-dark .bp3-running-text h2{
      color:#f5f8fa; }
  .bp3-running-text h3{
    color:#182026;
    font-weight:600;
    margin-top:40px;
    margin-bottom:20px; }
    .bp3-dark .bp3-running-text h3{
      color:#f5f8fa; }
  .bp3-running-text h4{
    color:#182026;
    font-weight:600;
    margin-top:40px;
    margin-bottom:20px; }
    .bp3-dark .bp3-running-text h4{
      color:#f5f8fa; }
  .bp3-running-text h5{
    color:#182026;
    font-weight:600;
    margin-top:40px;
    margin-bottom:20px; }
    .bp3-dark .bp3-running-text h5{
      color:#f5f8fa; }
  .bp3-running-text h6{
    color:#182026;
    font-weight:600;
    margin-top:40px;
    margin-bottom:20px; }
    .bp3-dark .bp3-running-text h6{
      color:#f5f8fa; }
  .bp3-running-text hr{
    margin:20px 0;
    border:none;
    border-bottom:1px solid rgba(16, 22, 26, 0.15); }
    .bp3-dark .bp3-running-text hr{
      border-color:rgba(255, 255, 255, 0.15); }
  .bp3-running-text p{
    margin:0 0 10px;
    padding:0; }

.bp3-text-large{
  font-size:16px; }

.bp3-text-small{
  font-size:12px; }
a{
  text-decoration:none;
  color:#106ba3; }
  a:hover{
    cursor:pointer;
    text-decoration:underline;
    color:#106ba3; }
  a .bp3-icon, a .bp3-icon-standard, a .bp3-icon-large{
    color:inherit; }
  a code,
  .bp3-dark a code{
    color:inherit; }
  .bp3-dark a,
  .bp3-dark a:hover{
    color:#48aff0; }
    .bp3-dark a .bp3-icon, .bp3-dark a .bp3-icon-standard, .bp3-dark a .bp3-icon-large,
    .bp3-dark a:hover .bp3-icon,
    .bp3-dark a:hover .bp3-icon-standard,
    .bp3-dark a:hover .bp3-icon-large{
      color:inherit; }
.bp3-running-text code, .bp3-code{
  text-transform:none;
  font-family:monospace;
  border-radius:3px;
  -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2);
          box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2);
  background:rgba(255, 255, 255, 0.7);
  padding:2px 5px;
  color:#5c7080;
  font-size:smaller; }
  .bp3-dark .bp3-running-text code, .bp3-running-text .bp3-dark code, .bp3-dark .bp3-code{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
    background:rgba(16, 22, 26, 0.3);
    color:#a7b6c2; }
  .bp3-running-text a > code, a > .bp3-code{
    color:#137cbd; }
    .bp3-dark .bp3-running-text a > code, .bp3-running-text .bp3-dark a > code, .bp3-dark a > .bp3-code{
      color:inherit; }

.bp3-running-text pre, .bp3-code-block{
  text-transform:none;
  font-family:monospace;
  display:block;
  margin:10px 0;
  border-radius:3px;
  -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.15);
          box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.15);
  background:rgba(255, 255, 255, 0.7);
  padding:13px 15px 12px;
  line-height:1.4;
  color:#182026;
  font-size:13px;
  word-break:break-all;
  word-wrap:break-word; }
  .bp3-dark .bp3-running-text pre, .bp3-running-text .bp3-dark pre, .bp3-dark .bp3-code-block{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
    background:rgba(16, 22, 26, 0.3);
    color:#f5f8fa; }
  .bp3-running-text pre > code, .bp3-code-block > code{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:none;
    padding:0;
    color:inherit;
    font-size:inherit; }

.bp3-running-text kbd, .bp3-key{
  display:-webkit-inline-box;
  display:-ms-inline-flexbox;
  display:inline-flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:center;
      -ms-flex-pack:center;
          justify-content:center;
  border-radius:3px;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
  background:#ffffff;
  min-width:24px;
  height:24px;
  padding:3px 6px;
  vertical-align:middle;
  line-height:24px;
  color:#5c7080;
  font-family:inherit;
  font-size:12px; }
  .bp3-running-text kbd .bp3-icon, .bp3-key .bp3-icon, .bp3-running-text kbd .bp3-icon-standard, .bp3-key .bp3-icon-standard, .bp3-running-text kbd .bp3-icon-large, .bp3-key .bp3-icon-large{
    margin-right:5px; }
  .bp3-dark .bp3-running-text kbd, .bp3-running-text .bp3-dark kbd, .bp3-dark .bp3-key{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4);
    background:#394b59;
    color:#a7b6c2; }
.bp3-running-text blockquote, .bp3-blockquote{
  margin:0 0 10px;
  border-left:solid 4px rgba(167, 182, 194, 0.5);
  padding:0 20px; }
  .bp3-dark .bp3-running-text blockquote, .bp3-running-text .bp3-dark blockquote, .bp3-dark .bp3-blockquote{
    border-color:rgba(115, 134, 148, 0.5); }
.bp3-running-text ul,
.bp3-running-text ol, .bp3-list{
  margin:10px 0;
  padding-left:30px; }
  .bp3-running-text ul li:not(:last-child), .bp3-running-text ol li:not(:last-child), .bp3-list li:not(:last-child){
    margin-bottom:5px; }
  .bp3-running-text ul ol, .bp3-running-text ol ol, .bp3-list ol,
  .bp3-running-text ul ul,
  .bp3-running-text ol ul,
  .bp3-list ul{
    margin-top:5px; }

.bp3-list-unstyled{
  margin:0;
  padding:0;
  list-style:none; }
  .bp3-list-unstyled li{
    padding:0; }
.bp3-rtl{
  text-align:right; }

.bp3-dark{
  color:#f5f8fa; }

:focus{
  outline:rgba(19, 124, 189, 0.6) auto 2px;
  outline-offset:2px;
  -moz-outline-radius:6px; }

.bp3-focus-disabled :focus{
  outline:none !important; }
  .bp3-focus-disabled :focus ~ .bp3-control-indicator{
    outline:none !important; }

.bp3-alert{
  max-width:400px;
  padding:20px; }

.bp3-alert-body{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex; }
  .bp3-alert-body .bp3-icon{
    margin-top:0;
    margin-right:20px;
    font-size:40px; }

.bp3-alert-footer{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:reverse;
      -ms-flex-direction:row-reverse;
          flex-direction:row-reverse;
  margin-top:10px; }
  .bp3-alert-footer .bp3-button{
    margin-left:10px; }
.bp3-breadcrumbs{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -ms-flex-wrap:wrap;
      flex-wrap:wrap;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  margin:0;
  cursor:default;
  height:30px;
  padding:0;
  list-style:none; }
  .bp3-breadcrumbs > li{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    -webkit-box-align:center;
        -ms-flex-align:center;
            align-items:center; }
    .bp3-breadcrumbs > li::after{
      display:block;
      margin:0 5px;
      background:url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16'%3e%3cpath fill-rule='evenodd' clip-rule='evenodd' d='M10.71 7.29l-4-4a1.003 1.003 0 0 0-1.42 1.42L8.59 8 5.3 11.29c-.19.18-.3.43-.3.71a1.003 1.003 0 0 0 1.71.71l4-4c.18-.18.29-.43.29-.71 0-.28-.11-.53-.29-.71z' fill='%235C7080'/%3e%3c/svg%3e");
      width:16px;
      height:16px;
      content:""; }
    .bp3-breadcrumbs > li:last-of-type::after{
      display:none; }

.bp3-breadcrumb,
.bp3-breadcrumb-current,
.bp3-breadcrumbs-collapsed{
  display:-webkit-inline-box;
  display:-ms-inline-flexbox;
  display:inline-flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  font-size:16px; }

.bp3-breadcrumb,
.bp3-breadcrumbs-collapsed{
  color:#5c7080; }

.bp3-breadcrumb:hover{
  text-decoration:none; }

.bp3-breadcrumb.bp3-disabled{
  cursor:not-allowed;
  color:rgba(92, 112, 128, 0.6); }

.bp3-breadcrumb .bp3-icon{
  margin-right:5px; }

.bp3-breadcrumb-current{
  color:inherit;
  font-weight:600; }
  .bp3-breadcrumb-current .bp3-input{
    vertical-align:baseline;
    font-size:inherit;
    font-weight:inherit; }

.bp3-breadcrumbs-collapsed{
  margin-right:2px;
  border:none;
  border-radius:3px;
  background:#ced9e0;
  cursor:pointer;
  padding:1px 5px;
  vertical-align:text-bottom; }
  .bp3-breadcrumbs-collapsed::before{
    display:block;
    background:url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16'%3e%3cg fill='%235C7080'%3e%3ccircle cx='2' cy='8.03' r='2'/%3e%3ccircle cx='14' cy='8.03' r='2'/%3e%3ccircle cx='8' cy='8.03' r='2'/%3e%3c/g%3e%3c/svg%3e") center no-repeat;
    width:16px;
    height:16px;
    content:""; }
  .bp3-breadcrumbs-collapsed:hover{
    background:#bfccd6;
    text-decoration:none;
    color:#182026; }

.bp3-dark .bp3-breadcrumb,
.bp3-dark .bp3-breadcrumbs-collapsed{
  color:#a7b6c2; }

.bp3-dark .bp3-breadcrumbs > li::after{
  color:#a7b6c2; }

.bp3-dark .bp3-breadcrumb.bp3-disabled{
  color:rgba(167, 182, 194, 0.6); }

.bp3-dark .bp3-breadcrumb-current{
  color:#f5f8fa; }

.bp3-dark .bp3-breadcrumbs-collapsed{
  background:rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-breadcrumbs-collapsed:hover{
    background:rgba(16, 22, 26, 0.6);
    color:#f5f8fa; }
.bp3-button{
  display:-webkit-inline-box;
  display:-ms-inline-flexbox;
  display:inline-flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:center;
      -ms-flex-pack:center;
          justify-content:center;
  border:none;
  border-radius:3px;
  cursor:pointer;
  padding:5px 10px;
  vertical-align:middle;
  text-align:left;
  font-size:14px;
  min-width:30px;
  min-height:30px; }
  .bp3-button > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-button > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-button::before,
  .bp3-button > *{
    margin-right:7px; }
  .bp3-button:empty::before,
  .bp3-button > :last-child{
    margin-right:0; }
  .bp3-button:empty{
    padding:0 !important; }
  .bp3-button:disabled, .bp3-button.bp3-disabled{
    cursor:not-allowed; }
  .bp3-button.bp3-fill{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    width:100%; }
  .bp3-button.bp3-align-right,
  .bp3-align-right .bp3-button{
    text-align:right; }
  .bp3-button.bp3-align-left,
  .bp3-align-left .bp3-button{
    text-align:left; }
  .bp3-button:not([class*="bp3-intent-"]){
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-color:#f5f8fa;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.8)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.8), rgba(255, 255, 255, 0));
    color:#182026; }
    .bp3-button:not([class*="bp3-intent-"]):hover{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
      background-clip:padding-box;
      background-color:#ebf1f5; }
    .bp3-button:not([class*="bp3-intent-"]):active, .bp3-button:not([class*="bp3-intent-"]).bp3-active{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#d8e1e8;
      background-image:none; }
    .bp3-button:not([class*="bp3-intent-"]):disabled, .bp3-button:not([class*="bp3-intent-"]).bp3-disabled{
      outline:none;
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(206, 217, 224, 0.5);
      background-image:none;
      cursor:not-allowed;
      color:rgba(92, 112, 128, 0.6); }
      .bp3-button:not([class*="bp3-intent-"]):disabled.bp3-active, .bp3-button:not([class*="bp3-intent-"]):disabled.bp3-active:hover, .bp3-button:not([class*="bp3-intent-"]).bp3-disabled.bp3-active, .bp3-button:not([class*="bp3-intent-"]).bp3-disabled.bp3-active:hover{
        background:rgba(206, 217, 224, 0.7); }
  .bp3-button.bp3-intent-primary{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#137cbd;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.1)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.1), rgba(255, 255, 255, 0));
    color:#ffffff; }
    .bp3-button.bp3-intent-primary:hover, .bp3-button.bp3-intent-primary:active, .bp3-button.bp3-intent-primary.bp3-active{
      color:#ffffff; }
    .bp3-button.bp3-intent-primary:hover{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
      background-color:#106ba3; }
    .bp3-button.bp3-intent-primary:active, .bp3-button.bp3-intent-primary.bp3-active{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#0e5a8a;
      background-image:none; }
    .bp3-button.bp3-intent-primary:disabled, .bp3-button.bp3-intent-primary.bp3-disabled{
      border-color:transparent;
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(19, 124, 189, 0.5);
      background-image:none;
      color:rgba(255, 255, 255, 0.6); }
  .bp3-button.bp3-intent-success{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#0f9960;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.1)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.1), rgba(255, 255, 255, 0));
    color:#ffffff; }
    .bp3-button.bp3-intent-success:hover, .bp3-button.bp3-intent-success:active, .bp3-button.bp3-intent-success.bp3-active{
      color:#ffffff; }
    .bp3-button.bp3-intent-success:hover{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
      background-color:#0d8050; }
    .bp3-button.bp3-intent-success:active, .bp3-button.bp3-intent-success.bp3-active{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#0a6640;
      background-image:none; }
    .bp3-button.bp3-intent-success:disabled, .bp3-button.bp3-intent-success.bp3-disabled{
      border-color:transparent;
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(15, 153, 96, 0.5);
      background-image:none;
      color:rgba(255, 255, 255, 0.6); }
  .bp3-button.bp3-intent-warning{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#d9822b;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.1)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.1), rgba(255, 255, 255, 0));
    color:#ffffff; }
    .bp3-button.bp3-intent-warning:hover, .bp3-button.bp3-intent-warning:active, .bp3-button.bp3-intent-warning.bp3-active{
      color:#ffffff; }
    .bp3-button.bp3-intent-warning:hover{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
      background-color:#bf7326; }
    .bp3-button.bp3-intent-warning:active, .bp3-button.bp3-intent-warning.bp3-active{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#a66321;
      background-image:none; }
    .bp3-button.bp3-intent-warning:disabled, .bp3-button.bp3-intent-warning.bp3-disabled{
      border-color:transparent;
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(217, 130, 43, 0.5);
      background-image:none;
      color:rgba(255, 255, 255, 0.6); }
  .bp3-button.bp3-intent-danger{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#db3737;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.1)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.1), rgba(255, 255, 255, 0));
    color:#ffffff; }
    .bp3-button.bp3-intent-danger:hover, .bp3-button.bp3-intent-danger:active, .bp3-button.bp3-intent-danger.bp3-active{
      color:#ffffff; }
    .bp3-button.bp3-intent-danger:hover{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
      background-color:#c23030; }
    .bp3-button.bp3-intent-danger:active, .bp3-button.bp3-intent-danger.bp3-active{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#a82a2a;
      background-image:none; }
    .bp3-button.bp3-intent-danger:disabled, .bp3-button.bp3-intent-danger.bp3-disabled{
      border-color:transparent;
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(219, 55, 55, 0.5);
      background-image:none;
      color:rgba(255, 255, 255, 0.6); }
  .bp3-button[class*="bp3-intent-"] .bp3-button-spinner .bp3-spinner-head{
    stroke:#ffffff; }
  .bp3-button.bp3-large,
  .bp3-large .bp3-button{
    min-width:40px;
    min-height:40px;
    padding:5px 15px;
    font-size:16px; }
    .bp3-button.bp3-large::before,
    .bp3-button.bp3-large > *,
    .bp3-large .bp3-button::before,
    .bp3-large .bp3-button > *{
      margin-right:10px; }
    .bp3-button.bp3-large:empty::before,
    .bp3-button.bp3-large > :last-child,
    .bp3-large .bp3-button:empty::before,
    .bp3-large .bp3-button > :last-child{
      margin-right:0; }
  .bp3-button.bp3-small,
  .bp3-small .bp3-button{
    min-width:24px;
    min-height:24px;
    padding:0 7px; }
  .bp3-button.bp3-loading{
    position:relative; }
    .bp3-button.bp3-loading[class*="bp3-icon-"]::before{
      visibility:hidden; }
    .bp3-button.bp3-loading .bp3-button-spinner{
      position:absolute;
      margin:0; }
    .bp3-button.bp3-loading > :not(.bp3-button-spinner){
      visibility:hidden; }
  .bp3-button[class*="bp3-icon-"]::before{
    line-height:1;
    font-family:"Icons16", sans-serif;
    font-size:16px;
    font-weight:400;
    font-style:normal;
    -moz-osx-font-smoothing:grayscale;
    -webkit-font-smoothing:antialiased;
    color:#5c7080; }
  .bp3-button .bp3-icon, .bp3-button .bp3-icon-standard, .bp3-button .bp3-icon-large{
    color:#5c7080; }
    .bp3-button .bp3-icon.bp3-align-right, .bp3-button .bp3-icon-standard.bp3-align-right, .bp3-button .bp3-icon-large.bp3-align-right{
      margin-left:7px; }
  .bp3-button .bp3-icon:first-child:last-child,
  .bp3-button .bp3-spinner + .bp3-icon:last-child{
    margin:0 -7px; }
  .bp3-dark .bp3-button:not([class*="bp3-intent-"]){
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
    background-color:#394b59;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.05)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.05), rgba(255, 255, 255, 0));
    color:#f5f8fa; }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"]):hover, .bp3-dark .bp3-button:not([class*="bp3-intent-"]):active, .bp3-dark .bp3-button:not([class*="bp3-intent-"]).bp3-active{
      color:#f5f8fa; }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"]):hover{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
      background-color:#30404d; }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"]):active, .bp3-dark .bp3-button:not([class*="bp3-intent-"]).bp3-active{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#202b33;
      background-image:none; }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"]):disabled, .bp3-dark .bp3-button:not([class*="bp3-intent-"]).bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(57, 75, 89, 0.5);
      background-image:none;
      color:rgba(167, 182, 194, 0.6); }
      .bp3-dark .bp3-button:not([class*="bp3-intent-"]):disabled.bp3-active, .bp3-dark .bp3-button:not([class*="bp3-intent-"]).bp3-disabled.bp3-active{
        background:rgba(57, 75, 89, 0.7); }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"]) .bp3-button-spinner .bp3-spinner-head{
      background:rgba(16, 22, 26, 0.5);
      stroke:#8a9ba8; }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"])[class*="bp3-icon-"]::before{
      color:#a7b6c2; }
    .bp3-dark .bp3-button:not([class*="bp3-intent-"]) .bp3-icon, .bp3-dark .bp3-button:not([class*="bp3-intent-"]) .bp3-icon-standard, .bp3-dark .bp3-button:not([class*="bp3-intent-"]) .bp3-icon-large{
      color:#a7b6c2; }
  .bp3-dark .bp3-button[class*="bp3-intent-"]{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-button[class*="bp3-intent-"]:hover{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-button[class*="bp3-intent-"]:active, .bp3-dark .bp3-button[class*="bp3-intent-"].bp3-active{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2); }
    .bp3-dark .bp3-button[class*="bp3-intent-"]:disabled, .bp3-dark .bp3-button[class*="bp3-intent-"].bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none;
      background-image:none;
      color:rgba(255, 255, 255, 0.3); }
    .bp3-dark .bp3-button[class*="bp3-intent-"] .bp3-button-spinner .bp3-spinner-head{
      stroke:#8a9ba8; }
  .bp3-button:disabled::before,
  .bp3-button:disabled .bp3-icon, .bp3-button:disabled .bp3-icon-standard, .bp3-button:disabled .bp3-icon-large, .bp3-button.bp3-disabled::before,
  .bp3-button.bp3-disabled .bp3-icon, .bp3-button.bp3-disabled .bp3-icon-standard, .bp3-button.bp3-disabled .bp3-icon-large, .bp3-button[class*="bp3-intent-"]::before,
  .bp3-button[class*="bp3-intent-"] .bp3-icon, .bp3-button[class*="bp3-intent-"] .bp3-icon-standard, .bp3-button[class*="bp3-intent-"] .bp3-icon-large{
    color:inherit !important; }
  .bp3-button.bp3-minimal{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:none; }
    .bp3-button.bp3-minimal:hover{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(167, 182, 194, 0.3);
      text-decoration:none;
      color:#182026; }
    .bp3-button.bp3-minimal:active, .bp3-button.bp3-minimal.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(115, 134, 148, 0.3);
      color:#182026; }
    .bp3-button.bp3-minimal:disabled, .bp3-button.bp3-minimal:disabled:hover, .bp3-button.bp3-minimal.bp3-disabled, .bp3-button.bp3-minimal.bp3-disabled:hover{
      background:none;
      cursor:not-allowed;
      color:rgba(92, 112, 128, 0.6); }
      .bp3-button.bp3-minimal:disabled.bp3-active, .bp3-button.bp3-minimal:disabled:hover.bp3-active, .bp3-button.bp3-minimal.bp3-disabled.bp3-active, .bp3-button.bp3-minimal.bp3-disabled:hover.bp3-active{
        background:rgba(115, 134, 148, 0.3); }
    .bp3-dark .bp3-button.bp3-minimal{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none;
      color:inherit; }
      .bp3-dark .bp3-button.bp3-minimal:hover, .bp3-dark .bp3-button.bp3-minimal:active, .bp3-dark .bp3-button.bp3-minimal.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none; }
      .bp3-dark .bp3-button.bp3-minimal:hover{
        background:rgba(138, 155, 168, 0.15); }
      .bp3-dark .bp3-button.bp3-minimal:active, .bp3-dark .bp3-button.bp3-minimal.bp3-active{
        background:rgba(138, 155, 168, 0.3);
        color:#f5f8fa; }
      .bp3-dark .bp3-button.bp3-minimal:disabled, .bp3-dark .bp3-button.bp3-minimal:disabled:hover, .bp3-dark .bp3-button.bp3-minimal.bp3-disabled, .bp3-dark .bp3-button.bp3-minimal.bp3-disabled:hover{
        background:none;
        cursor:not-allowed;
        color:rgba(167, 182, 194, 0.6); }
        .bp3-dark .bp3-button.bp3-minimal:disabled.bp3-active, .bp3-dark .bp3-button.bp3-minimal:disabled:hover.bp3-active, .bp3-dark .bp3-button.bp3-minimal.bp3-disabled.bp3-active, .bp3-dark .bp3-button.bp3-minimal.bp3-disabled:hover.bp3-active{
          background:rgba(138, 155, 168, 0.3); }
    .bp3-button.bp3-minimal.bp3-intent-primary{
      color:#106ba3; }
      .bp3-button.bp3-minimal.bp3-intent-primary:hover, .bp3-button.bp3-minimal.bp3-intent-primary:active, .bp3-button.bp3-minimal.bp3-intent-primary.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#106ba3; }
      .bp3-button.bp3-minimal.bp3-intent-primary:hover{
        background:rgba(19, 124, 189, 0.15);
        color:#106ba3; }
      .bp3-button.bp3-minimal.bp3-intent-primary:active, .bp3-button.bp3-minimal.bp3-intent-primary.bp3-active{
        background:rgba(19, 124, 189, 0.3);
        color:#106ba3; }
      .bp3-button.bp3-minimal.bp3-intent-primary:disabled, .bp3-button.bp3-minimal.bp3-intent-primary.bp3-disabled{
        background:none;
        color:rgba(16, 107, 163, 0.5); }
        .bp3-button.bp3-minimal.bp3-intent-primary:disabled.bp3-active, .bp3-button.bp3-minimal.bp3-intent-primary.bp3-disabled.bp3-active{
          background:rgba(19, 124, 189, 0.3); }
      .bp3-button.bp3-minimal.bp3-intent-primary .bp3-button-spinner .bp3-spinner-head{
        stroke:#106ba3; }
      .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary{
        color:#48aff0; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary:hover{
          background:rgba(19, 124, 189, 0.2);
          color:#48aff0; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary:active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary.bp3-active{
          background:rgba(19, 124, 189, 0.3);
          color:#48aff0; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary:disabled, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary.bp3-disabled{
          background:none;
          color:rgba(72, 175, 240, 0.5); }
          .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary:disabled.bp3-active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-primary.bp3-disabled.bp3-active{
            background:rgba(19, 124, 189, 0.3); }
    .bp3-button.bp3-minimal.bp3-intent-success{
      color:#0d8050; }
      .bp3-button.bp3-minimal.bp3-intent-success:hover, .bp3-button.bp3-minimal.bp3-intent-success:active, .bp3-button.bp3-minimal.bp3-intent-success.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#0d8050; }
      .bp3-button.bp3-minimal.bp3-intent-success:hover{
        background:rgba(15, 153, 96, 0.15);
        color:#0d8050; }
      .bp3-button.bp3-minimal.bp3-intent-success:active, .bp3-button.bp3-minimal.bp3-intent-success.bp3-active{
        background:rgba(15, 153, 96, 0.3);
        color:#0d8050; }
      .bp3-button.bp3-minimal.bp3-intent-success:disabled, .bp3-button.bp3-minimal.bp3-intent-success.bp3-disabled{
        background:none;
        color:rgba(13, 128, 80, 0.5); }
        .bp3-button.bp3-minimal.bp3-intent-success:disabled.bp3-active, .bp3-button.bp3-minimal.bp3-intent-success.bp3-disabled.bp3-active{
          background:rgba(15, 153, 96, 0.3); }
      .bp3-button.bp3-minimal.bp3-intent-success .bp3-button-spinner .bp3-spinner-head{
        stroke:#0d8050; }
      .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success{
        color:#3dcc91; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success:hover{
          background:rgba(15, 153, 96, 0.2);
          color:#3dcc91; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success:active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success.bp3-active{
          background:rgba(15, 153, 96, 0.3);
          color:#3dcc91; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success:disabled, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success.bp3-disabled{
          background:none;
          color:rgba(61, 204, 145, 0.5); }
          .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success:disabled.bp3-active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-success.bp3-disabled.bp3-active{
            background:rgba(15, 153, 96, 0.3); }
    .bp3-button.bp3-minimal.bp3-intent-warning{
      color:#bf7326; }
      .bp3-button.bp3-minimal.bp3-intent-warning:hover, .bp3-button.bp3-minimal.bp3-intent-warning:active, .bp3-button.bp3-minimal.bp3-intent-warning.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#bf7326; }
      .bp3-button.bp3-minimal.bp3-intent-warning:hover{
        background:rgba(217, 130, 43, 0.15);
        color:#bf7326; }
      .bp3-button.bp3-minimal.bp3-intent-warning:active, .bp3-button.bp3-minimal.bp3-intent-warning.bp3-active{
        background:rgba(217, 130, 43, 0.3);
        color:#bf7326; }
      .bp3-button.bp3-minimal.bp3-intent-warning:disabled, .bp3-button.bp3-minimal.bp3-intent-warning.bp3-disabled{
        background:none;
        color:rgba(191, 115, 38, 0.5); }
        .bp3-button.bp3-minimal.bp3-intent-warning:disabled.bp3-active, .bp3-button.bp3-minimal.bp3-intent-warning.bp3-disabled.bp3-active{
          background:rgba(217, 130, 43, 0.3); }
      .bp3-button.bp3-minimal.bp3-intent-warning .bp3-button-spinner .bp3-spinner-head{
        stroke:#bf7326; }
      .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning{
        color:#ffb366; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning:hover{
          background:rgba(217, 130, 43, 0.2);
          color:#ffb366; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning:active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning.bp3-active{
          background:rgba(217, 130, 43, 0.3);
          color:#ffb366; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning:disabled, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning.bp3-disabled{
          background:none;
          color:rgba(255, 179, 102, 0.5); }
          .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning:disabled.bp3-active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-warning.bp3-disabled.bp3-active{
            background:rgba(217, 130, 43, 0.3); }
    .bp3-button.bp3-minimal.bp3-intent-danger{
      color:#c23030; }
      .bp3-button.bp3-minimal.bp3-intent-danger:hover, .bp3-button.bp3-minimal.bp3-intent-danger:active, .bp3-button.bp3-minimal.bp3-intent-danger.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#c23030; }
      .bp3-button.bp3-minimal.bp3-intent-danger:hover{
        background:rgba(219, 55, 55, 0.15);
        color:#c23030; }
      .bp3-button.bp3-minimal.bp3-intent-danger:active, .bp3-button.bp3-minimal.bp3-intent-danger.bp3-active{
        background:rgba(219, 55, 55, 0.3);
        color:#c23030; }
      .bp3-button.bp3-minimal.bp3-intent-danger:disabled, .bp3-button.bp3-minimal.bp3-intent-danger.bp3-disabled{
        background:none;
        color:rgba(194, 48, 48, 0.5); }
        .bp3-button.bp3-minimal.bp3-intent-danger:disabled.bp3-active, .bp3-button.bp3-minimal.bp3-intent-danger.bp3-disabled.bp3-active{
          background:rgba(219, 55, 55, 0.3); }
      .bp3-button.bp3-minimal.bp3-intent-danger .bp3-button-spinner .bp3-spinner-head{
        stroke:#c23030; }
      .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger{
        color:#ff7373; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger:hover{
          background:rgba(219, 55, 55, 0.2);
          color:#ff7373; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger:active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger.bp3-active{
          background:rgba(219, 55, 55, 0.3);
          color:#ff7373; }
        .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger:disabled, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger.bp3-disabled{
          background:none;
          color:rgba(255, 115, 115, 0.5); }
          .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger:disabled.bp3-active, .bp3-dark .bp3-button.bp3-minimal.bp3-intent-danger.bp3-disabled.bp3-active{
            background:rgba(219, 55, 55, 0.3); }

a.bp3-button{
  text-align:center;
  text-decoration:none;
  -webkit-transition:none;
  transition:none; }
  a.bp3-button, a.bp3-button:hover, a.bp3-button:active{
    color:#182026; }
  a.bp3-button.bp3-disabled{
    color:rgba(92, 112, 128, 0.6); }

.bp3-button-text{
  -webkit-box-flex:0;
      -ms-flex:0 1 auto;
          flex:0 1 auto; }

.bp3-button.bp3-align-left .bp3-button-text, .bp3-button.bp3-align-right .bp3-button-text,
.bp3-button-group.bp3-align-left .bp3-button-text,
.bp3-button-group.bp3-align-right .bp3-button-text{
  -webkit-box-flex:1;
      -ms-flex:1 1 auto;
          flex:1 1 auto; }
.bp3-button-group{
  display:-webkit-inline-box;
  display:-ms-inline-flexbox;
  display:inline-flex; }
  .bp3-button-group .bp3-button{
    -webkit-box-flex:0;
        -ms-flex:0 0 auto;
            flex:0 0 auto;
    position:relative;
    z-index:4; }
    .bp3-button-group .bp3-button:focus{
      z-index:5; }
    .bp3-button-group .bp3-button:hover{
      z-index:6; }
    .bp3-button-group .bp3-button:active, .bp3-button-group .bp3-button.bp3-active{
      z-index:7; }
    .bp3-button-group .bp3-button:disabled, .bp3-button-group .bp3-button.bp3-disabled{
      z-index:3; }
    .bp3-button-group .bp3-button[class*="bp3-intent-"]{
      z-index:9; }
      .bp3-button-group .bp3-button[class*="bp3-intent-"]:focus{
        z-index:10; }
      .bp3-button-group .bp3-button[class*="bp3-intent-"]:hover{
        z-index:11; }
      .bp3-button-group .bp3-button[class*="bp3-intent-"]:active, .bp3-button-group .bp3-button[class*="bp3-intent-"].bp3-active{
        z-index:12; }
      .bp3-button-group .bp3-button[class*="bp3-intent-"]:disabled, .bp3-button-group .bp3-button[class*="bp3-intent-"].bp3-disabled{
        z-index:8; }
  .bp3-button-group:not(.bp3-minimal) > .bp3-popover-wrapper:not(:first-child) .bp3-button,
  .bp3-button-group:not(.bp3-minimal) > .bp3-button:not(:first-child){
    border-top-left-radius:0;
    border-bottom-left-radius:0; }
  .bp3-button-group:not(.bp3-minimal) > .bp3-popover-wrapper:not(:last-child) .bp3-button,
  .bp3-button-group:not(.bp3-minimal) > .bp3-button:not(:last-child){
    margin-right:-1px;
    border-top-right-radius:0;
    border-bottom-right-radius:0; }
  .bp3-button-group.bp3-minimal .bp3-button{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:none; }
    .bp3-button-group.bp3-minimal .bp3-button:hover{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(167, 182, 194, 0.3);
      text-decoration:none;
      color:#182026; }
    .bp3-button-group.bp3-minimal .bp3-button:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(115, 134, 148, 0.3);
      color:#182026; }
    .bp3-button-group.bp3-minimal .bp3-button:disabled, .bp3-button-group.bp3-minimal .bp3-button:disabled:hover, .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled, .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled:hover{
      background:none;
      cursor:not-allowed;
      color:rgba(92, 112, 128, 0.6); }
      .bp3-button-group.bp3-minimal .bp3-button:disabled.bp3-active, .bp3-button-group.bp3-minimal .bp3-button:disabled:hover.bp3-active, .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled.bp3-active, .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled:hover.bp3-active{
        background:rgba(115, 134, 148, 0.3); }
    .bp3-dark .bp3-button-group.bp3-minimal .bp3-button{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none;
      color:inherit; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:hover, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:hover{
        background:rgba(138, 155, 168, 0.15); }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-active{
        background:rgba(138, 155, 168, 0.3);
        color:#f5f8fa; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:disabled, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:disabled:hover, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled:hover{
        background:none;
        cursor:not-allowed;
        color:rgba(167, 182, 194, 0.6); }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:disabled.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button:disabled:hover.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-disabled:hover.bp3-active{
          background:rgba(138, 155, 168, 0.3); }
    .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary{
      color:#106ba3; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:hover, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#106ba3; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:hover{
        background:rgba(19, 124, 189, 0.15);
        color:#106ba3; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-active{
        background:rgba(19, 124, 189, 0.3);
        color:#106ba3; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:disabled, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-disabled{
        background:none;
        color:rgba(16, 107, 163, 0.5); }
        .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:disabled.bp3-active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-disabled.bp3-active{
          background:rgba(19, 124, 189, 0.3); }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary .bp3-button-spinner .bp3-spinner-head{
        stroke:#106ba3; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary{
        color:#48aff0; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:hover{
          background:rgba(19, 124, 189, 0.2);
          color:#48aff0; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-active{
          background:rgba(19, 124, 189, 0.3);
          color:#48aff0; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:disabled, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-disabled{
          background:none;
          color:rgba(72, 175, 240, 0.5); }
          .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary:disabled.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-primary.bp3-disabled.bp3-active{
            background:rgba(19, 124, 189, 0.3); }
    .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success{
      color:#0d8050; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:hover, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#0d8050; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:hover{
        background:rgba(15, 153, 96, 0.15);
        color:#0d8050; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-active{
        background:rgba(15, 153, 96, 0.3);
        color:#0d8050; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:disabled, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-disabled{
        background:none;
        color:rgba(13, 128, 80, 0.5); }
        .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:disabled.bp3-active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-disabled.bp3-active{
          background:rgba(15, 153, 96, 0.3); }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success .bp3-button-spinner .bp3-spinner-head{
        stroke:#0d8050; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success{
        color:#3dcc91; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:hover{
          background:rgba(15, 153, 96, 0.2);
          color:#3dcc91; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-active{
          background:rgba(15, 153, 96, 0.3);
          color:#3dcc91; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:disabled, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-disabled{
          background:none;
          color:rgba(61, 204, 145, 0.5); }
          .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success:disabled.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-success.bp3-disabled.bp3-active{
            background:rgba(15, 153, 96, 0.3); }
    .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning{
      color:#bf7326; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:hover, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#bf7326; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:hover{
        background:rgba(217, 130, 43, 0.15);
        color:#bf7326; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-active{
        background:rgba(217, 130, 43, 0.3);
        color:#bf7326; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:disabled, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-disabled{
        background:none;
        color:rgba(191, 115, 38, 0.5); }
        .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:disabled.bp3-active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-disabled.bp3-active{
          background:rgba(217, 130, 43, 0.3); }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning .bp3-button-spinner .bp3-spinner-head{
        stroke:#bf7326; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning{
        color:#ffb366; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:hover{
          background:rgba(217, 130, 43, 0.2);
          color:#ffb366; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-active{
          background:rgba(217, 130, 43, 0.3);
          color:#ffb366; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:disabled, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-disabled{
          background:none;
          color:rgba(255, 179, 102, 0.5); }
          .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning:disabled.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-warning.bp3-disabled.bp3-active{
            background:rgba(217, 130, 43, 0.3); }
    .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger{
      color:#c23030; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:hover, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-active{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:none;
        color:#c23030; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:hover{
        background:rgba(219, 55, 55, 0.15);
        color:#c23030; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-active{
        background:rgba(219, 55, 55, 0.3);
        color:#c23030; }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:disabled, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-disabled{
        background:none;
        color:rgba(194, 48, 48, 0.5); }
        .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:disabled.bp3-active, .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-disabled.bp3-active{
          background:rgba(219, 55, 55, 0.3); }
      .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger .bp3-button-spinner .bp3-spinner-head{
        stroke:#c23030; }
      .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger{
        color:#ff7373; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:hover{
          background:rgba(219, 55, 55, 0.2);
          color:#ff7373; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-active{
          background:rgba(219, 55, 55, 0.3);
          color:#ff7373; }
        .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:disabled, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-disabled{
          background:none;
          color:rgba(255, 115, 115, 0.5); }
          .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger:disabled.bp3-active, .bp3-dark .bp3-button-group.bp3-minimal .bp3-button.bp3-intent-danger.bp3-disabled.bp3-active{
            background:rgba(219, 55, 55, 0.3); }
  .bp3-button-group .bp3-popover-wrapper,
  .bp3-button-group .bp3-popover-target{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto; }
  .bp3-button-group.bp3-fill{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    width:100%; }
  .bp3-button-group .bp3-button.bp3-fill,
  .bp3-button-group.bp3-fill .bp3-button:not(.bp3-fixed){
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto; }
  .bp3-button-group.bp3-vertical{
    -webkit-box-orient:vertical;
    -webkit-box-direction:normal;
        -ms-flex-direction:column;
            flex-direction:column;
    -webkit-box-align:stretch;
        -ms-flex-align:stretch;
            align-items:stretch;
    vertical-align:top; }
    .bp3-button-group.bp3-vertical.bp3-fill{
      width:unset;
      height:100%; }
    .bp3-button-group.bp3-vertical .bp3-button{
      margin-right:0 !important;
      width:100%; }
    .bp3-button-group.bp3-vertical:not(.bp3-minimal) > .bp3-popover-wrapper:first-child .bp3-button,
    .bp3-button-group.bp3-vertical:not(.bp3-minimal) > .bp3-button:first-child{
      border-radius:3px 3px 0 0; }
    .bp3-button-group.bp3-vertical:not(.bp3-minimal) > .bp3-popover-wrapper:last-child .bp3-button,
    .bp3-button-group.bp3-vertical:not(.bp3-minimal) > .bp3-button:last-child{
      border-radius:0 0 3px 3px; }
    .bp3-button-group.bp3-vertical:not(.bp3-minimal) > .bp3-popover-wrapper:not(:last-child) .bp3-button,
    .bp3-button-group.bp3-vertical:not(.bp3-minimal) > .bp3-button:not(:last-child){
      margin-bottom:-1px; }
  .bp3-button-group.bp3-align-left .bp3-button{
    text-align:left; }
  .bp3-dark .bp3-button-group:not(.bp3-minimal) > .bp3-popover-wrapper:not(:last-child) .bp3-button,
  .bp3-dark .bp3-button-group:not(.bp3-minimal) > .bp3-button:not(:last-child){
    margin-right:1px; }
  .bp3-dark .bp3-button-group.bp3-vertical > .bp3-popover-wrapper:not(:last-child) .bp3-button,
  .bp3-dark .bp3-button-group.bp3-vertical > .bp3-button:not(:last-child){
    margin-bottom:1px; }
.bp3-callout{
  line-height:1.5;
  font-size:14px;
  position:relative;
  border-radius:3px;
  background-color:rgba(138, 155, 168, 0.15);
  width:100%;
  padding:10px 12px 9px; }
  .bp3-callout[class*="bp3-icon-"]{
    padding-left:40px; }
    .bp3-callout[class*="bp3-icon-"]::before{
      line-height:1;
      font-family:"Icons20", sans-serif;
      font-size:20px;
      font-weight:400;
      font-style:normal;
      -moz-osx-font-smoothing:grayscale;
      -webkit-font-smoothing:antialiased;
      position:absolute;
      top:10px;
      left:10px;
      color:#5c7080; }
  .bp3-callout.bp3-callout-icon{
    padding-left:40px; }
    .bp3-callout.bp3-callout-icon > .bp3-icon:first-child{
      position:absolute;
      top:10px;
      left:10px;
      color:#5c7080; }
  .bp3-callout .bp3-heading{
    margin-top:0;
    margin-bottom:5px;
    line-height:20px; }
    .bp3-callout .bp3-heading:last-child{
      margin-bottom:0; }
  .bp3-dark .bp3-callout{
    background-color:rgba(138, 155, 168, 0.2); }
    .bp3-dark .bp3-callout[class*="bp3-icon-"]::before{
      color:#a7b6c2; }
  .bp3-callout.bp3-intent-primary{
    background-color:rgba(19, 124, 189, 0.15); }
    .bp3-callout.bp3-intent-primary[class*="bp3-icon-"]::before,
    .bp3-callout.bp3-intent-primary > .bp3-icon:first-child,
    .bp3-callout.bp3-intent-primary .bp3-heading{
      color:#106ba3; }
    .bp3-dark .bp3-callout.bp3-intent-primary{
      background-color:rgba(19, 124, 189, 0.25); }
      .bp3-dark .bp3-callout.bp3-intent-primary[class*="bp3-icon-"]::before,
      .bp3-dark .bp3-callout.bp3-intent-primary > .bp3-icon:first-child,
      .bp3-dark .bp3-callout.bp3-intent-primary .bp3-heading{
        color:#48aff0; }
  .bp3-callout.bp3-intent-success{
    background-color:rgba(15, 153, 96, 0.15); }
    .bp3-callout.bp3-intent-success[class*="bp3-icon-"]::before,
    .bp3-callout.bp3-intent-success > .bp3-icon:first-child,
    .bp3-callout.bp3-intent-success .bp3-heading{
      color:#0d8050; }
    .bp3-dark .bp3-callout.bp3-intent-success{
      background-color:rgba(15, 153, 96, 0.25); }
      .bp3-dark .bp3-callout.bp3-intent-success[class*="bp3-icon-"]::before,
      .bp3-dark .bp3-callout.bp3-intent-success > .bp3-icon:first-child,
      .bp3-dark .bp3-callout.bp3-intent-success .bp3-heading{
        color:#3dcc91; }
  .bp3-callout.bp3-intent-warning{
    background-color:rgba(217, 130, 43, 0.15); }
    .bp3-callout.bp3-intent-warning[class*="bp3-icon-"]::before,
    .bp3-callout.bp3-intent-warning > .bp3-icon:first-child,
    .bp3-callout.bp3-intent-warning .bp3-heading{
      color:#bf7326; }
    .bp3-dark .bp3-callout.bp3-intent-warning{
      background-color:rgba(217, 130, 43, 0.25); }
      .bp3-dark .bp3-callout.bp3-intent-warning[class*="bp3-icon-"]::before,
      .bp3-dark .bp3-callout.bp3-intent-warning > .bp3-icon:first-child,
      .bp3-dark .bp3-callout.bp3-intent-warning .bp3-heading{
        color:#ffb366; }
  .bp3-callout.bp3-intent-danger{
    background-color:rgba(219, 55, 55, 0.15); }
    .bp3-callout.bp3-intent-danger[class*="bp3-icon-"]::before,
    .bp3-callout.bp3-intent-danger > .bp3-icon:first-child,
    .bp3-callout.bp3-intent-danger .bp3-heading{
      color:#c23030; }
    .bp3-dark .bp3-callout.bp3-intent-danger{
      background-color:rgba(219, 55, 55, 0.25); }
      .bp3-dark .bp3-callout.bp3-intent-danger[class*="bp3-icon-"]::before,
      .bp3-dark .bp3-callout.bp3-intent-danger > .bp3-icon:first-child,
      .bp3-dark .bp3-callout.bp3-intent-danger .bp3-heading{
        color:#ff7373; }
  .bp3-running-text .bp3-callout{
    margin:20px 0; }
.bp3-card{
  border-radius:3px;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.15), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.15), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0);
  background-color:#ffffff;
  padding:20px;
  -webkit-transition:-webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:-webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), box-shadow 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), box-shadow 200ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 200ms cubic-bezier(0.4, 1, 0.75, 0.9); }
  .bp3-card.bp3-dark,
  .bp3-dark .bp3-card{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0);
    background-color:#30404d; }

.bp3-elevation-0{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.15), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.15), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0); }
  .bp3-elevation-0.bp3-dark,
  .bp3-dark .bp3-elevation-0{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), 0 0 0 rgba(16, 22, 26, 0), 0 0 0 rgba(16, 22, 26, 0); }

.bp3-elevation-1{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-elevation-1.bp3-dark,
  .bp3-dark .bp3-elevation-1{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4); }

.bp3-elevation-2{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 1px 1px rgba(16, 22, 26, 0.2), 0 2px 6px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 1px 1px rgba(16, 22, 26, 0.2), 0 2px 6px rgba(16, 22, 26, 0.2); }
  .bp3-elevation-2.bp3-dark,
  .bp3-dark .bp3-elevation-2{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.4), 0 2px 6px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.4), 0 2px 6px rgba(16, 22, 26, 0.4); }

.bp3-elevation-3{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2); }
  .bp3-elevation-3.bp3-dark,
  .bp3-dark .bp3-elevation-3{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4); }

.bp3-elevation-4{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2); }
  .bp3-elevation-4.bp3-dark,
  .bp3-dark .bp3-elevation-4{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4); }

.bp3-card.bp3-interactive:hover{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
  cursor:pointer; }
  .bp3-card.bp3-interactive:hover.bp3-dark,
  .bp3-dark .bp3-card.bp3-interactive:hover{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4); }

.bp3-card.bp3-interactive:active{
  opacity:0.9;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
  -webkit-transition-duration:0;
          transition-duration:0; }
  .bp3-card.bp3-interactive:active.bp3-dark,
  .bp3-dark .bp3-card.bp3-interactive:active{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4); }

.bp3-collapse{
  height:0;
  overflow-y:hidden;
  -webkit-transition:height 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:height 200ms cubic-bezier(0.4, 1, 0.75, 0.9); }
  .bp3-collapse .bp3-collapse-body{
    -webkit-transition:-webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:-webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9); }
    .bp3-collapse .bp3-collapse-body[aria-hidden="true"]{
      display:none; }

.bp3-context-menu .bp3-popover-target{
  display:block; }

.bp3-context-menu-popover-target{
  position:fixed; }

.bp3-divider{
  margin:5px;
  border-right:1px solid rgba(16, 22, 26, 0.15);
  border-bottom:1px solid rgba(16, 22, 26, 0.15); }
  .bp3-dark .bp3-divider{
    border-color:rgba(16, 22, 26, 0.4); }
.bp3-dialog-container{
  opacity:1;
  -webkit-transform:scale(1);
          transform:scale(1);
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:center;
      -ms-flex-pack:center;
          justify-content:center;
  width:100%;
  min-height:100%;
  pointer-events:none;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-dialog-container.bp3-overlay-enter > .bp3-dialog, .bp3-dialog-container.bp3-overlay-appear > .bp3-dialog{
    opacity:0;
    -webkit-transform:scale(0.5);
            transform:scale(0.5); }
  .bp3-dialog-container.bp3-overlay-enter-active > .bp3-dialog, .bp3-dialog-container.bp3-overlay-appear-active > .bp3-dialog{
    opacity:1;
    -webkit-transform:scale(1);
            transform:scale(1);
    -webkit-transition-property:opacity, -webkit-transform;
    transition-property:opacity, -webkit-transform;
    transition-property:opacity, transform;
    transition-property:opacity, transform, -webkit-transform;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
            transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-dialog-container.bp3-overlay-exit > .bp3-dialog{
    opacity:1;
    -webkit-transform:scale(1);
            transform:scale(1); }
  .bp3-dialog-container.bp3-overlay-exit-active > .bp3-dialog{
    opacity:0;
    -webkit-transform:scale(0.5);
            transform:scale(0.5);
    -webkit-transition-property:opacity, -webkit-transform;
    transition-property:opacity, -webkit-transform;
    transition-property:opacity, transform;
    transition-property:opacity, transform, -webkit-transform;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
            transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
    -webkit-transition-delay:0;
            transition-delay:0; }

.bp3-dialog{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:vertical;
  -webkit-box-direction:normal;
      -ms-flex-direction:column;
          flex-direction:column;
  margin:30px 0;
  border-radius:6px;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
  background:#ebf1f5;
  width:500px;
  padding-bottom:20px;
  pointer-events:all;
  -webkit-user-select:text;
     -moz-user-select:text;
      -ms-user-select:text;
          user-select:text; }
  .bp3-dialog:focus{
    outline:0; }
  .bp3-dialog.bp3-dark,
  .bp3-dark .bp3-dialog{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
    background:#293742;
    color:#f5f8fa; }

.bp3-dialog-header{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  border-radius:6px 6px 0 0;
  -webkit-box-shadow:0 1px 0 rgba(16, 22, 26, 0.15);
          box-shadow:0 1px 0 rgba(16, 22, 26, 0.15);
  background:#ffffff;
  min-height:40px;
  padding-right:5px;
  padding-left:20px; }
  .bp3-dialog-header .bp3-icon-large,
  .bp3-dialog-header .bp3-icon{
    -webkit-box-flex:0;
        -ms-flex:0 0 auto;
            flex:0 0 auto;
    margin-right:10px;
    color:#5c7080; }
  .bp3-dialog-header .bp3-heading{
    overflow:hidden;
    text-overflow:ellipsis;
    white-space:nowrap;
    word-wrap:normal;
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto;
    margin:0;
    line-height:inherit; }
    .bp3-dialog-header .bp3-heading:last-child{
      margin-right:20px; }
  .bp3-dark .bp3-dialog-header{
    -webkit-box-shadow:0 1px 0 rgba(16, 22, 26, 0.4);
            box-shadow:0 1px 0 rgba(16, 22, 26, 0.4);
    background:#30404d; }
    .bp3-dark .bp3-dialog-header .bp3-icon-large,
    .bp3-dark .bp3-dialog-header .bp3-icon{
      color:#a7b6c2; }

.bp3-dialog-body{
  -webkit-box-flex:1;
      -ms-flex:1 1 auto;
          flex:1 1 auto;
  margin:20px;
  line-height:18px; }

.bp3-dialog-footer{
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  margin:0 20px; }

.bp3-dialog-footer-actions{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-pack:end;
      -ms-flex-pack:end;
          justify-content:flex-end; }
  .bp3-dialog-footer-actions .bp3-button{
    margin-left:10px; }
.bp3-drawer{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:vertical;
  -webkit-box-direction:normal;
      -ms-flex-direction:column;
          flex-direction:column;
  margin:0;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
  background:#ffffff;
  padding:0; }
  .bp3-drawer:focus{
    outline:0; }
  .bp3-drawer.bp3-position-top{
    top:0;
    right:0;
    left:0;
    height:50%; }
    .bp3-drawer.bp3-position-top.bp3-overlay-enter, .bp3-drawer.bp3-position-top.bp3-overlay-appear{
      -webkit-transform:translateY(-100%);
              transform:translateY(-100%); }
    .bp3-drawer.bp3-position-top.bp3-overlay-enter-active, .bp3-drawer.bp3-position-top.bp3-overlay-appear-active{
      -webkit-transform:translateY(0);
              transform:translateY(0);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:200ms;
              transition-duration:200ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
    .bp3-drawer.bp3-position-top.bp3-overlay-exit{
      -webkit-transform:translateY(0);
              transform:translateY(0); }
    .bp3-drawer.bp3-position-top.bp3-overlay-exit-active{
      -webkit-transform:translateY(-100%);
              transform:translateY(-100%);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:100ms;
              transition-duration:100ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
  .bp3-drawer.bp3-position-bottom{
    right:0;
    bottom:0;
    left:0;
    height:50%; }
    .bp3-drawer.bp3-position-bottom.bp3-overlay-enter, .bp3-drawer.bp3-position-bottom.bp3-overlay-appear{
      -webkit-transform:translateY(100%);
              transform:translateY(100%); }
    .bp3-drawer.bp3-position-bottom.bp3-overlay-enter-active, .bp3-drawer.bp3-position-bottom.bp3-overlay-appear-active{
      -webkit-transform:translateY(0);
              transform:translateY(0);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:200ms;
              transition-duration:200ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
    .bp3-drawer.bp3-position-bottom.bp3-overlay-exit{
      -webkit-transform:translateY(0);
              transform:translateY(0); }
    .bp3-drawer.bp3-position-bottom.bp3-overlay-exit-active{
      -webkit-transform:translateY(100%);
              transform:translateY(100%);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:100ms;
              transition-duration:100ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
  .bp3-drawer.bp3-position-left{
    top:0;
    bottom:0;
    left:0;
    width:50%; }
    .bp3-drawer.bp3-position-left.bp3-overlay-enter, .bp3-drawer.bp3-position-left.bp3-overlay-appear{
      -webkit-transform:translateX(-100%);
              transform:translateX(-100%); }
    .bp3-drawer.bp3-position-left.bp3-overlay-enter-active, .bp3-drawer.bp3-position-left.bp3-overlay-appear-active{
      -webkit-transform:translateX(0);
              transform:translateX(0);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:200ms;
              transition-duration:200ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
    .bp3-drawer.bp3-position-left.bp3-overlay-exit{
      -webkit-transform:translateX(0);
              transform:translateX(0); }
    .bp3-drawer.bp3-position-left.bp3-overlay-exit-active{
      -webkit-transform:translateX(-100%);
              transform:translateX(-100%);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:100ms;
              transition-duration:100ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
  .bp3-drawer.bp3-position-right{
    top:0;
    right:0;
    bottom:0;
    width:50%; }
    .bp3-drawer.bp3-position-right.bp3-overlay-enter, .bp3-drawer.bp3-position-right.bp3-overlay-appear{
      -webkit-transform:translateX(100%);
              transform:translateX(100%); }
    .bp3-drawer.bp3-position-right.bp3-overlay-enter-active, .bp3-drawer.bp3-position-right.bp3-overlay-appear-active{
      -webkit-transform:translateX(0);
              transform:translateX(0);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:200ms;
              transition-duration:200ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
    .bp3-drawer.bp3-position-right.bp3-overlay-exit{
      -webkit-transform:translateX(0);
              transform:translateX(0); }
    .bp3-drawer.bp3-position-right.bp3-overlay-exit-active{
      -webkit-transform:translateX(100%);
              transform:translateX(100%);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:100ms;
              transition-duration:100ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
  .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
  .bp3-position-right):not(.bp3-vertical){
    top:0;
    right:0;
    bottom:0;
    width:50%; }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right):not(.bp3-vertical).bp3-overlay-enter, .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right):not(.bp3-vertical).bp3-overlay-appear{
      -webkit-transform:translateX(100%);
              transform:translateX(100%); }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right):not(.bp3-vertical).bp3-overlay-enter-active, .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right):not(.bp3-vertical).bp3-overlay-appear-active{
      -webkit-transform:translateX(0);
              transform:translateX(0);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:200ms;
              transition-duration:200ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right):not(.bp3-vertical).bp3-overlay-exit{
      -webkit-transform:translateX(0);
              transform:translateX(0); }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right):not(.bp3-vertical).bp3-overlay-exit-active{
      -webkit-transform:translateX(100%);
              transform:translateX(100%);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:100ms;
              transition-duration:100ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
  .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
  .bp3-position-right).bp3-vertical{
    right:0;
    bottom:0;
    left:0;
    height:50%; }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right).bp3-vertical.bp3-overlay-enter, .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right).bp3-vertical.bp3-overlay-appear{
      -webkit-transform:translateY(100%);
              transform:translateY(100%); }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right).bp3-vertical.bp3-overlay-enter-active, .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right).bp3-vertical.bp3-overlay-appear-active{
      -webkit-transform:translateY(0);
              transform:translateY(0);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:200ms;
              transition-duration:200ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right).bp3-vertical.bp3-overlay-exit{
      -webkit-transform:translateY(0);
              transform:translateY(0); }
    .bp3-drawer:not(.bp3-position-top):not(.bp3-position-bottom):not(.bp3-position-left):not(
    .bp3-position-right).bp3-vertical.bp3-overlay-exit-active{
      -webkit-transform:translateY(100%);
              transform:translateY(100%);
      -webkit-transition-property:-webkit-transform;
      transition-property:-webkit-transform;
      transition-property:transform;
      transition-property:transform, -webkit-transform;
      -webkit-transition-duration:100ms;
              transition-duration:100ms;
      -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
              transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
      -webkit-transition-delay:0;
              transition-delay:0; }
  .bp3-drawer.bp3-dark,
  .bp3-dark .bp3-drawer{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
    background:#30404d;
    color:#f5f8fa; }

.bp3-drawer-header{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  position:relative;
  border-radius:0;
  -webkit-box-shadow:0 1px 0 rgba(16, 22, 26, 0.15);
          box-shadow:0 1px 0 rgba(16, 22, 26, 0.15);
  min-height:40px;
  padding:5px;
  padding-left:20px; }
  .bp3-drawer-header .bp3-icon-large,
  .bp3-drawer-header .bp3-icon{
    -webkit-box-flex:0;
        -ms-flex:0 0 auto;
            flex:0 0 auto;
    margin-right:10px;
    color:#5c7080; }
  .bp3-drawer-header .bp3-heading{
    overflow:hidden;
    text-overflow:ellipsis;
    white-space:nowrap;
    word-wrap:normal;
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto;
    margin:0;
    line-height:inherit; }
    .bp3-drawer-header .bp3-heading:last-child{
      margin-right:20px; }
  .bp3-dark .bp3-drawer-header{
    -webkit-box-shadow:0 1px 0 rgba(16, 22, 26, 0.4);
            box-shadow:0 1px 0 rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-drawer-header .bp3-icon-large,
    .bp3-dark .bp3-drawer-header .bp3-icon{
      color:#a7b6c2; }

.bp3-drawer-body{
  -webkit-box-flex:1;
      -ms-flex:1 1 auto;
          flex:1 1 auto;
  overflow:auto;
  line-height:18px; }

.bp3-drawer-footer{
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  position:relative;
  -webkit-box-shadow:inset 0 1px 0 rgba(16, 22, 26, 0.15);
          box-shadow:inset 0 1px 0 rgba(16, 22, 26, 0.15);
  padding:10px 20px; }
  .bp3-dark .bp3-drawer-footer{
    -webkit-box-shadow:inset 0 1px 0 rgba(16, 22, 26, 0.4);
            box-shadow:inset 0 1px 0 rgba(16, 22, 26, 0.4); }
.bp3-editable-text{
  display:inline-block;
  position:relative;
  cursor:text;
  max-width:100%;
  vertical-align:top;
  white-space:nowrap; }
  .bp3-editable-text::before{
    position:absolute;
    top:-3px;
    right:-3px;
    bottom:-3px;
    left:-3px;
    border-radius:3px;
    content:"";
    -webkit-transition:background-color 100ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:background-color 100ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:background-color 100ms cubic-bezier(0.4, 1, 0.75, 0.9), box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:background-color 100ms cubic-bezier(0.4, 1, 0.75, 0.9), box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9); }
  .bp3-editable-text:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.15);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.15); }
  .bp3-editable-text.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
    background-color:#ffffff; }
  .bp3-editable-text.bp3-disabled::before{
    -webkit-box-shadow:none;
            box-shadow:none; }
  .bp3-editable-text.bp3-intent-primary .bp3-editable-text-input,
  .bp3-editable-text.bp3-intent-primary .bp3-editable-text-content{
    color:#137cbd; }
  .bp3-editable-text.bp3-intent-primary:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(19, 124, 189, 0.4);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(19, 124, 189, 0.4); }
  .bp3-editable-text.bp3-intent-primary.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-editable-text.bp3-intent-success .bp3-editable-text-input,
  .bp3-editable-text.bp3-intent-success .bp3-editable-text-content{
    color:#0f9960; }
  .bp3-editable-text.bp3-intent-success:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px rgba(15, 153, 96, 0.4);
            box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px rgba(15, 153, 96, 0.4); }
  .bp3-editable-text.bp3-intent-success.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-editable-text.bp3-intent-warning .bp3-editable-text-input,
  .bp3-editable-text.bp3-intent-warning .bp3-editable-text-content{
    color:#d9822b; }
  .bp3-editable-text.bp3-intent-warning:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px rgba(217, 130, 43, 0.4);
            box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px rgba(217, 130, 43, 0.4); }
  .bp3-editable-text.bp3-intent-warning.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-editable-text.bp3-intent-danger .bp3-editable-text-input,
  .bp3-editable-text.bp3-intent-danger .bp3-editable-text-content{
    color:#db3737; }
  .bp3-editable-text.bp3-intent-danger:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px rgba(219, 55, 55, 0.4);
            box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px rgba(219, 55, 55, 0.4); }
  .bp3-editable-text.bp3-intent-danger.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-dark .bp3-editable-text:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(255, 255, 255, 0.15);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(255, 255, 255, 0.15); }
  .bp3-dark .bp3-editable-text.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
    background-color:rgba(16, 22, 26, 0.3); }
  .bp3-dark .bp3-editable-text.bp3-disabled::before{
    -webkit-box-shadow:none;
            box-shadow:none; }
  .bp3-dark .bp3-editable-text.bp3-intent-primary .bp3-editable-text-content{
    color:#48aff0; }
  .bp3-dark .bp3-editable-text.bp3-intent-primary:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(72, 175, 240, 0), 0 0 0 0 rgba(72, 175, 240, 0), inset 0 0 0 1px rgba(72, 175, 240, 0.4);
            box-shadow:0 0 0 0 rgba(72, 175, 240, 0), 0 0 0 0 rgba(72, 175, 240, 0), inset 0 0 0 1px rgba(72, 175, 240, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-primary.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #48aff0, 0 0 0 3px rgba(72, 175, 240, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px #48aff0, 0 0 0 3px rgba(72, 175, 240, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-success .bp3-editable-text-content{
    color:#3dcc91; }
  .bp3-dark .bp3-editable-text.bp3-intent-success:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(61, 204, 145, 0), 0 0 0 0 rgba(61, 204, 145, 0), inset 0 0 0 1px rgba(61, 204, 145, 0.4);
            box-shadow:0 0 0 0 rgba(61, 204, 145, 0), 0 0 0 0 rgba(61, 204, 145, 0), inset 0 0 0 1px rgba(61, 204, 145, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-success.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #3dcc91, 0 0 0 3px rgba(61, 204, 145, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px #3dcc91, 0 0 0 3px rgba(61, 204, 145, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-warning .bp3-editable-text-content{
    color:#ffb366; }
  .bp3-dark .bp3-editable-text.bp3-intent-warning:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(255, 179, 102, 0), 0 0 0 0 rgba(255, 179, 102, 0), inset 0 0 0 1px rgba(255, 179, 102, 0.4);
            box-shadow:0 0 0 0 rgba(255, 179, 102, 0), 0 0 0 0 rgba(255, 179, 102, 0), inset 0 0 0 1px rgba(255, 179, 102, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-warning.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #ffb366, 0 0 0 3px rgba(255, 179, 102, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px #ffb366, 0 0 0 3px rgba(255, 179, 102, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-danger .bp3-editable-text-content{
    color:#ff7373; }
  .bp3-dark .bp3-editable-text.bp3-intent-danger:hover::before{
    -webkit-box-shadow:0 0 0 0 rgba(255, 115, 115, 0), 0 0 0 0 rgba(255, 115, 115, 0), inset 0 0 0 1px rgba(255, 115, 115, 0.4);
            box-shadow:0 0 0 0 rgba(255, 115, 115, 0), 0 0 0 0 rgba(255, 115, 115, 0), inset 0 0 0 1px rgba(255, 115, 115, 0.4); }
  .bp3-dark .bp3-editable-text.bp3-intent-danger.bp3-editable-text-editing::before{
    -webkit-box-shadow:0 0 0 1px #ff7373, 0 0 0 3px rgba(255, 115, 115, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px #ff7373, 0 0 0 3px rgba(255, 115, 115, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }

.bp3-editable-text-input,
.bp3-editable-text-content{
  display:inherit;
  position:relative;
  min-width:inherit;
  max-width:inherit;
  vertical-align:top;
  text-transform:inherit;
  letter-spacing:inherit;
  color:inherit;
  font:inherit;
  resize:none; }

.bp3-editable-text-input{
  border:none;
  -webkit-box-shadow:none;
          box-shadow:none;
  background:none;
  width:100%;
  padding:0;
  white-space:pre-wrap; }
  .bp3-editable-text-input::-webkit-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-editable-text-input::-moz-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-editable-text-input:-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-editable-text-input::-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-editable-text-input::placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-editable-text-input:focus{
    outline:none; }
  .bp3-editable-text-input::-ms-clear{
    display:none; }

.bp3-editable-text-content{
  overflow:hidden;
  padding-right:2px;
  text-overflow:ellipsis;
  white-space:pre; }
  .bp3-editable-text-editing > .bp3-editable-text-content{
    position:absolute;
    left:0;
    visibility:hidden; }
  .bp3-editable-text-placeholder > .bp3-editable-text-content{
    color:rgba(92, 112, 128, 0.6); }
    .bp3-dark .bp3-editable-text-placeholder > .bp3-editable-text-content{
      color:rgba(167, 182, 194, 0.6); }

.bp3-editable-text.bp3-multiline{
  display:block; }
  .bp3-editable-text.bp3-multiline .bp3-editable-text-content{
    overflow:auto;
    white-space:pre-wrap;
    word-wrap:break-word; }
.bp3-control-group{
  -webkit-transform:translateZ(0);
          transform:translateZ(0);
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:stretch;
      -ms-flex-align:stretch;
          align-items:stretch; }
  .bp3-control-group > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-control-group > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-control-group .bp3-button,
  .bp3-control-group .bp3-html-select,
  .bp3-control-group .bp3-input,
  .bp3-control-group .bp3-select{
    position:relative; }
  .bp3-control-group .bp3-input{
    z-index:2;
    border-radius:inherit; }
    .bp3-control-group .bp3-input:focus{
      z-index:14;
      border-radius:3px; }
    .bp3-control-group .bp3-input[class*="bp3-intent"]{
      z-index:13; }
      .bp3-control-group .bp3-input[class*="bp3-intent"]:focus{
        z-index:15; }
    .bp3-control-group .bp3-input[readonly], .bp3-control-group .bp3-input:disabled, .bp3-control-group .bp3-input.bp3-disabled{
      z-index:1; }
  .bp3-control-group .bp3-input-group[class*="bp3-intent"] .bp3-input{
    z-index:13; }
    .bp3-control-group .bp3-input-group[class*="bp3-intent"] .bp3-input:focus{
      z-index:15; }
  .bp3-control-group .bp3-button,
  .bp3-control-group .bp3-html-select select,
  .bp3-control-group .bp3-select select{
    -webkit-transform:translateZ(0);
            transform:translateZ(0);
    z-index:4;
    border-radius:inherit; }
    .bp3-control-group .bp3-button:focus,
    .bp3-control-group .bp3-html-select select:focus,
    .bp3-control-group .bp3-select select:focus{
      z-index:5; }
    .bp3-control-group .bp3-button:hover,
    .bp3-control-group .bp3-html-select select:hover,
    .bp3-control-group .bp3-select select:hover{
      z-index:6; }
    .bp3-control-group .bp3-button:active,
    .bp3-control-group .bp3-html-select select:active,
    .bp3-control-group .bp3-select select:active{
      z-index:7; }
    .bp3-control-group .bp3-button[readonly], .bp3-control-group .bp3-button:disabled, .bp3-control-group .bp3-button.bp3-disabled,
    .bp3-control-group .bp3-html-select select[readonly],
    .bp3-control-group .bp3-html-select select:disabled,
    .bp3-control-group .bp3-html-select select.bp3-disabled,
    .bp3-control-group .bp3-select select[readonly],
    .bp3-control-group .bp3-select select:disabled,
    .bp3-control-group .bp3-select select.bp3-disabled{
      z-index:3; }
    .bp3-control-group .bp3-button[class*="bp3-intent"],
    .bp3-control-group .bp3-html-select select[class*="bp3-intent"],
    .bp3-control-group .bp3-select select[class*="bp3-intent"]{
      z-index:9; }
      .bp3-control-group .bp3-button[class*="bp3-intent"]:focus,
      .bp3-control-group .bp3-html-select select[class*="bp3-intent"]:focus,
      .bp3-control-group .bp3-select select[class*="bp3-intent"]:focus{
        z-index:10; }
      .bp3-control-group .bp3-button[class*="bp3-intent"]:hover,
      .bp3-control-group .bp3-html-select select[class*="bp3-intent"]:hover,
      .bp3-control-group .bp3-select select[class*="bp3-intent"]:hover{
        z-index:11; }
      .bp3-control-group .bp3-button[class*="bp3-intent"]:active,
      .bp3-control-group .bp3-html-select select[class*="bp3-intent"]:active,
      .bp3-control-group .bp3-select select[class*="bp3-intent"]:active{
        z-index:12; }
      .bp3-control-group .bp3-button[class*="bp3-intent"][readonly], .bp3-control-group .bp3-button[class*="bp3-intent"]:disabled, .bp3-control-group .bp3-button[class*="bp3-intent"].bp3-disabled,
      .bp3-control-group .bp3-html-select select[class*="bp3-intent"][readonly],
      .bp3-control-group .bp3-html-select select[class*="bp3-intent"]:disabled,
      .bp3-control-group .bp3-html-select select[class*="bp3-intent"].bp3-disabled,
      .bp3-control-group .bp3-select select[class*="bp3-intent"][readonly],
      .bp3-control-group .bp3-select select[class*="bp3-intent"]:disabled,
      .bp3-control-group .bp3-select select[class*="bp3-intent"].bp3-disabled{
        z-index:8; }
  .bp3-control-group .bp3-input-group > .bp3-icon,
  .bp3-control-group .bp3-input-group > .bp3-button,
  .bp3-control-group .bp3-input-group > .bp3-input-action{
    z-index:16; }
  .bp3-control-group .bp3-select::after,
  .bp3-control-group .bp3-html-select::after,
  .bp3-control-group .bp3-select > .bp3-icon,
  .bp3-control-group .bp3-html-select > .bp3-icon{
    z-index:17; }
  .bp3-control-group:not(.bp3-vertical) > *{
    margin-right:-1px; }
  .bp3-dark .bp3-control-group:not(.bp3-vertical) > *{
    margin-right:0; }
  .bp3-dark .bp3-control-group:not(.bp3-vertical) > .bp3-button + .bp3-button{
    margin-left:1px; }
  .bp3-control-group .bp3-popover-wrapper,
  .bp3-control-group .bp3-popover-target{
    border-radius:inherit; }
  .bp3-control-group > :first-child{
    border-radius:3px 0 0 3px; }
  .bp3-control-group > :last-child{
    margin-right:0;
    border-radius:0 3px 3px 0; }
  .bp3-control-group > :only-child{
    margin-right:0;
    border-radius:3px; }
  .bp3-control-group .bp3-input-group .bp3-button{
    border-radius:3px; }
  .bp3-control-group > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto; }
  .bp3-control-group.bp3-fill > *:not(.bp3-fixed){
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto; }
  .bp3-control-group.bp3-vertical{
    -webkit-box-orient:vertical;
    -webkit-box-direction:normal;
        -ms-flex-direction:column;
            flex-direction:column; }
    .bp3-control-group.bp3-vertical > *{
      margin-top:-1px; }
    .bp3-control-group.bp3-vertical > :first-child{
      margin-top:0;
      border-radius:3px 3px 0 0; }
    .bp3-control-group.bp3-vertical > :last-child{
      border-radius:0 0 3px 3px; }
.bp3-control{
  display:block;
  position:relative;
  margin-bottom:10px;
  cursor:pointer;
  text-transform:none; }
  .bp3-control input:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#137cbd;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.1)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.1), rgba(255, 255, 255, 0));
    color:#ffffff; }
  .bp3-control:hover input:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#106ba3; }
  .bp3-control input:not(:disabled):active:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background:#0e5a8a; }
  .bp3-control input:disabled:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(19, 124, 189, 0.5); }
  .bp3-dark .bp3-control input:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-control:hover input:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
    background-color:#106ba3; }
  .bp3-dark .bp3-control input:not(:disabled):active:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#0e5a8a; }
  .bp3-dark .bp3-control input:disabled:checked ~ .bp3-control-indicator{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(14, 90, 138, 0.5); }
  .bp3-control:not(.bp3-align-right){
    padding-left:26px; }
    .bp3-control:not(.bp3-align-right) .bp3-control-indicator{
      margin-left:-26px; }
  .bp3-control.bp3-align-right{
    padding-right:26px; }
    .bp3-control.bp3-align-right .bp3-control-indicator{
      margin-right:-26px; }
  .bp3-control.bp3-disabled{
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-control.bp3-inline{
    display:inline-block;
    margin-right:20px; }
  .bp3-control input{
    position:absolute;
    top:0;
    left:0;
    opacity:0;
    z-index:-1; }
  .bp3-control .bp3-control-indicator{
    display:inline-block;
    position:relative;
    margin-top:-3px;
    margin-right:10px;
    border:none;
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-clip:padding-box;
    background-color:#f5f8fa;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.8)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.8), rgba(255, 255, 255, 0));
    cursor:pointer;
    width:1em;
    height:1em;
    vertical-align:middle;
    font-size:16px;
    -webkit-user-select:none;
       -moz-user-select:none;
        -ms-user-select:none;
            user-select:none; }
    .bp3-control .bp3-control-indicator::before{
      display:block;
      width:1em;
      height:1em;
      content:""; }
  .bp3-control:hover .bp3-control-indicator{
    background-color:#ebf1f5; }
  .bp3-control input:not(:disabled):active ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background:#d8e1e8; }
  .bp3-control input:disabled ~ .bp3-control-indicator{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(206, 217, 224, 0.5);
    cursor:not-allowed; }
  .bp3-control input:focus ~ .bp3-control-indicator{
    outline:rgba(19, 124, 189, 0.6) auto 2px;
    outline-offset:2px;
    -moz-outline-radius:6px; }
  .bp3-control.bp3-align-right .bp3-control-indicator{
    float:right;
    margin-top:1px;
    margin-left:10px; }
  .bp3-control.bp3-large{
    font-size:16px; }
    .bp3-control.bp3-large:not(.bp3-align-right){
      padding-left:30px; }
      .bp3-control.bp3-large:not(.bp3-align-right) .bp3-control-indicator{
        margin-left:-30px; }
    .bp3-control.bp3-large.bp3-align-right{
      padding-right:30px; }
      .bp3-control.bp3-large.bp3-align-right .bp3-control-indicator{
        margin-right:-30px; }
    .bp3-control.bp3-large .bp3-control-indicator{
      font-size:20px; }
    .bp3-control.bp3-large.bp3-align-right .bp3-control-indicator{
      margin-top:0; }
  .bp3-control.bp3-checkbox input:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#137cbd;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.1)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.1), rgba(255, 255, 255, 0));
    color:#ffffff; }
  .bp3-control.bp3-checkbox:hover input:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 -1px 0 rgba(16, 22, 26, 0.2);
    background-color:#106ba3; }
  .bp3-control.bp3-checkbox input:not(:disabled):active:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background:#0e5a8a; }
  .bp3-control.bp3-checkbox input:disabled:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(19, 124, 189, 0.5); }
  .bp3-dark .bp3-control.bp3-checkbox input:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-control.bp3-checkbox:hover input:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
    background-color:#106ba3; }
  .bp3-dark .bp3-control.bp3-checkbox input:not(:disabled):active:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#0e5a8a; }
  .bp3-dark .bp3-control.bp3-checkbox input:disabled:indeterminate ~ .bp3-control-indicator{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(14, 90, 138, 0.5); }
  .bp3-control.bp3-checkbox .bp3-control-indicator{
    border-radius:3px; }
  .bp3-control.bp3-checkbox input:checked ~ .bp3-control-indicator::before{
    background-image:url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16'%3e%3cpath fill-rule='evenodd' clip-rule='evenodd' d='M12 5c-.28 0-.53.11-.71.29L7 9.59l-2.29-2.3a1.003 1.003 0 0 0-1.42 1.42l3 3c.18.18.43.29.71.29s.53-.11.71-.29l5-5A1.003 1.003 0 0 0 12 5z' fill='white'/%3e%3c/svg%3e"); }
  .bp3-control.bp3-checkbox input:indeterminate ~ .bp3-control-indicator::before{
    background-image:url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16'%3e%3cpath fill-rule='evenodd' clip-rule='evenodd' d='M11 7H5c-.55 0-1 .45-1 1s.45 1 1 1h6c.55 0 1-.45 1-1s-.45-1-1-1z' fill='white'/%3e%3c/svg%3e"); }
  .bp3-control.bp3-radio .bp3-control-indicator{
    border-radius:50%; }
  .bp3-control.bp3-radio input:checked ~ .bp3-control-indicator::before{
    background-image:radial-gradient(#ffffff, #ffffff 28%, transparent 32%); }
  .bp3-control.bp3-radio input:checked:disabled ~ .bp3-control-indicator::before{
    opacity:0.5; }
  .bp3-control.bp3-radio input:focus ~ .bp3-control-indicator{
    -moz-outline-radius:16px; }
  .bp3-control.bp3-switch input ~ .bp3-control-indicator{
    background:rgba(167, 182, 194, 0.5); }
  .bp3-control.bp3-switch:hover input ~ .bp3-control-indicator{
    background:rgba(115, 134, 148, 0.5); }
  .bp3-control.bp3-switch input:not(:disabled):active ~ .bp3-control-indicator{
    background:rgba(92, 112, 128, 0.5); }
  .bp3-control.bp3-switch input:disabled ~ .bp3-control-indicator{
    background:rgba(206, 217, 224, 0.5); }
    .bp3-control.bp3-switch input:disabled ~ .bp3-control-indicator::before{
      background:rgba(255, 255, 255, 0.8); }
  .bp3-control.bp3-switch input:checked ~ .bp3-control-indicator{
    background:#137cbd; }
  .bp3-control.bp3-switch:hover input:checked ~ .bp3-control-indicator{
    background:#106ba3; }
  .bp3-control.bp3-switch input:checked:not(:disabled):active ~ .bp3-control-indicator{
    background:#0e5a8a; }
  .bp3-control.bp3-switch input:checked:disabled ~ .bp3-control-indicator{
    background:rgba(19, 124, 189, 0.5); }
    .bp3-control.bp3-switch input:checked:disabled ~ .bp3-control-indicator::before{
      background:rgba(255, 255, 255, 0.8); }
  .bp3-control.bp3-switch:not(.bp3-align-right){
    padding-left:38px; }
    .bp3-control.bp3-switch:not(.bp3-align-right) .bp3-control-indicator{
      margin-left:-38px; }
  .bp3-control.bp3-switch.bp3-align-right{
    padding-right:38px; }
    .bp3-control.bp3-switch.bp3-align-right .bp3-control-indicator{
      margin-right:-38px; }
  .bp3-control.bp3-switch .bp3-control-indicator{
    border:none;
    border-radius:1.75em;
    -webkit-box-shadow:none !important;
            box-shadow:none !important;
    width:auto;
    min-width:1.75em;
    -webkit-transition:background-color 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:background-color 100ms cubic-bezier(0.4, 1, 0.75, 0.9); }
    .bp3-control.bp3-switch .bp3-control-indicator::before{
      position:absolute;
      left:0;
      margin:2px;
      border-radius:50%;
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.2);
      background:#ffffff;
      width:calc(1em - 4px);
      height:calc(1em - 4px);
      -webkit-transition:left 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
      transition:left 100ms cubic-bezier(0.4, 1, 0.75, 0.9); }
  .bp3-control.bp3-switch input:checked ~ .bp3-control-indicator::before{
    left:calc(100% - 1em); }
  .bp3-control.bp3-switch.bp3-large:not(.bp3-align-right){
    padding-left:45px; }
    .bp3-control.bp3-switch.bp3-large:not(.bp3-align-right) .bp3-control-indicator{
      margin-left:-45px; }
  .bp3-control.bp3-switch.bp3-large.bp3-align-right{
    padding-right:45px; }
    .bp3-control.bp3-switch.bp3-large.bp3-align-right .bp3-control-indicator{
      margin-right:-45px; }
  .bp3-dark .bp3-control.bp3-switch input ~ .bp3-control-indicator{
    background:rgba(16, 22, 26, 0.5); }
  .bp3-dark .bp3-control.bp3-switch:hover input ~ .bp3-control-indicator{
    background:rgba(16, 22, 26, 0.7); }
  .bp3-dark .bp3-control.bp3-switch input:not(:disabled):active ~ .bp3-control-indicator{
    background:rgba(16, 22, 26, 0.9); }
  .bp3-dark .bp3-control.bp3-switch input:disabled ~ .bp3-control-indicator{
    background:rgba(57, 75, 89, 0.5); }
    .bp3-dark .bp3-control.bp3-switch input:disabled ~ .bp3-control-indicator::before{
      background:rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-control.bp3-switch input:checked ~ .bp3-control-indicator{
    background:#137cbd; }
  .bp3-dark .bp3-control.bp3-switch:hover input:checked ~ .bp3-control-indicator{
    background:#106ba3; }
  .bp3-dark .bp3-control.bp3-switch input:checked:not(:disabled):active ~ .bp3-control-indicator{
    background:#0e5a8a; }
  .bp3-dark .bp3-control.bp3-switch input:checked:disabled ~ .bp3-control-indicator{
    background:rgba(14, 90, 138, 0.5); }
    .bp3-dark .bp3-control.bp3-switch input:checked:disabled ~ .bp3-control-indicator::before{
      background:rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-control.bp3-switch .bp3-control-indicator::before{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
    background:#394b59; }
  .bp3-dark .bp3-control.bp3-switch input:checked ~ .bp3-control-indicator::before{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4); }
  .bp3-control.bp3-switch .bp3-switch-inner-text{
    text-align:center;
    font-size:0.7em; }
  .bp3-control.bp3-switch .bp3-control-indicator-child:first-child{
    visibility:hidden;
    margin-right:1.2em;
    margin-left:0.5em;
    line-height:0; }
  .bp3-control.bp3-switch .bp3-control-indicator-child:last-child{
    visibility:visible;
    margin-right:0.5em;
    margin-left:1.2em;
    line-height:1em; }
  .bp3-control.bp3-switch input:checked ~ .bp3-control-indicator .bp3-control-indicator-child:first-child{
    visibility:visible;
    line-height:1em; }
  .bp3-control.bp3-switch input:checked ~ .bp3-control-indicator .bp3-control-indicator-child:last-child{
    visibility:hidden;
    line-height:0; }
  .bp3-dark .bp3-control{
    color:#f5f8fa; }
    .bp3-dark .bp3-control.bp3-disabled{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-control .bp3-control-indicator{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
      background-color:#394b59;
      background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.05)), to(rgba(255, 255, 255, 0)));
      background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.05), rgba(255, 255, 255, 0)); }
    .bp3-dark .bp3-control:hover .bp3-control-indicator{
      background-color:#30404d; }
    .bp3-dark .bp3-control input:not(:disabled):active ~ .bp3-control-indicator{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background:#202b33; }
    .bp3-dark .bp3-control input:disabled ~ .bp3-control-indicator{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(57, 75, 89, 0.5);
      cursor:not-allowed; }
    .bp3-dark .bp3-control.bp3-checkbox input:disabled:checked ~ .bp3-control-indicator, .bp3-dark .bp3-control.bp3-checkbox input:disabled:indeterminate ~ .bp3-control-indicator{
      color:rgba(167, 182, 194, 0.6); }
.bp3-file-input{
  display:inline-block;
  position:relative;
  cursor:pointer;
  height:30px; }
  .bp3-file-input input{
    opacity:0;
    margin:0;
    min-width:200px; }
    .bp3-file-input input:disabled + .bp3-file-upload-input,
    .bp3-file-input input.bp3-disabled + .bp3-file-upload-input{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(206, 217, 224, 0.5);
      cursor:not-allowed;
      color:rgba(92, 112, 128, 0.6);
      resize:none; }
      .bp3-file-input input:disabled + .bp3-file-upload-input::after,
      .bp3-file-input input.bp3-disabled + .bp3-file-upload-input::after{
        outline:none;
        -webkit-box-shadow:none;
                box-shadow:none;
        background-color:rgba(206, 217, 224, 0.5);
        background-image:none;
        cursor:not-allowed;
        color:rgba(92, 112, 128, 0.6); }
        .bp3-file-input input:disabled + .bp3-file-upload-input::after.bp3-active, .bp3-file-input input:disabled + .bp3-file-upload-input::after.bp3-active:hover,
        .bp3-file-input input.bp3-disabled + .bp3-file-upload-input::after.bp3-active,
        .bp3-file-input input.bp3-disabled + .bp3-file-upload-input::after.bp3-active:hover{
          background:rgba(206, 217, 224, 0.7); }
      .bp3-dark .bp3-file-input input:disabled + .bp3-file-upload-input, .bp3-dark
      .bp3-file-input input.bp3-disabled + .bp3-file-upload-input{
        -webkit-box-shadow:none;
                box-shadow:none;
        background:rgba(57, 75, 89, 0.5);
        color:rgba(167, 182, 194, 0.6); }
        .bp3-dark .bp3-file-input input:disabled + .bp3-file-upload-input::after, .bp3-dark
        .bp3-file-input input.bp3-disabled + .bp3-file-upload-input::after{
          -webkit-box-shadow:none;
                  box-shadow:none;
          background-color:rgba(57, 75, 89, 0.5);
          background-image:none;
          color:rgba(167, 182, 194, 0.6); }
          .bp3-dark .bp3-file-input input:disabled + .bp3-file-upload-input::after.bp3-active, .bp3-dark
          .bp3-file-input input.bp3-disabled + .bp3-file-upload-input::after.bp3-active{
            background:rgba(57, 75, 89, 0.7); }
  .bp3-file-input.bp3-file-input-has-selection .bp3-file-upload-input{
    color:#182026; }
  .bp3-dark .bp3-file-input.bp3-file-input-has-selection .bp3-file-upload-input{
    color:#f5f8fa; }
  .bp3-file-input.bp3-fill{
    width:100%; }
  .bp3-file-input.bp3-large,
  .bp3-large .bp3-file-input{
    height:40px; }
  .bp3-file-input .bp3-file-upload-input-custom-text::after{
    content:attr(bp3-button-text); }

.bp3-file-upload-input{
  outline:none;
  border:none;
  border-radius:3px;
  -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
  background:#ffffff;
  height:30px;
  padding:0 10px;
  vertical-align:middle;
  line-height:30px;
  color:#182026;
  font-size:14px;
  font-weight:400;
  -webkit-transition:-webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:-webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  -webkit-appearance:none;
     -moz-appearance:none;
          appearance:none;
  overflow:hidden;
  text-overflow:ellipsis;
  white-space:nowrap;
  word-wrap:normal;
  position:absolute;
  top:0;
  right:0;
  left:0;
  padding-right:80px;
  color:rgba(92, 112, 128, 0.6);
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-file-upload-input::-webkit-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-file-upload-input::-moz-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-file-upload-input:-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-file-upload-input::-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-file-upload-input::placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-file-upload-input:focus, .bp3-file-upload-input.bp3-active{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-file-upload-input[type="search"], .bp3-file-upload-input.bp3-round{
    border-radius:30px;
    -webkit-box-sizing:border-box;
            box-sizing:border-box;
    padding-left:10px; }
  .bp3-file-upload-input[readonly]{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.15);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.15); }
  .bp3-file-upload-input:disabled, .bp3-file-upload-input.bp3-disabled{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(206, 217, 224, 0.5);
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6);
    resize:none; }
  .bp3-file-upload-input::after{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-color:#f5f8fa;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.8)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.8), rgba(255, 255, 255, 0));
    color:#182026;
    min-width:24px;
    min-height:24px;
    overflow:hidden;
    text-overflow:ellipsis;
    white-space:nowrap;
    word-wrap:normal;
    position:absolute;
    top:0;
    right:0;
    margin:3px;
    border-radius:3px;
    width:70px;
    text-align:center;
    line-height:24px;
    content:"Browse"; }
    .bp3-file-upload-input::after:hover{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
      background-clip:padding-box;
      background-color:#ebf1f5; }
    .bp3-file-upload-input::after:active, .bp3-file-upload-input::after.bp3-active{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#d8e1e8;
      background-image:none; }
    .bp3-file-upload-input::after:disabled, .bp3-file-upload-input::after.bp3-disabled{
      outline:none;
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(206, 217, 224, 0.5);
      background-image:none;
      cursor:not-allowed;
      color:rgba(92, 112, 128, 0.6); }
      .bp3-file-upload-input::after:disabled.bp3-active, .bp3-file-upload-input::after:disabled.bp3-active:hover, .bp3-file-upload-input::after.bp3-disabled.bp3-active, .bp3-file-upload-input::after.bp3-disabled.bp3-active:hover{
        background:rgba(206, 217, 224, 0.7); }
  .bp3-file-upload-input:hover::after{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-clip:padding-box;
    background-color:#ebf1f5; }
  .bp3-file-upload-input:active::after{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#d8e1e8;
    background-image:none; }
  .bp3-large .bp3-file-upload-input{
    height:40px;
    line-height:40px;
    font-size:16px;
    padding-right:95px; }
    .bp3-large .bp3-file-upload-input[type="search"], .bp3-large .bp3-file-upload-input.bp3-round{
      padding:0 15px; }
    .bp3-large .bp3-file-upload-input::after{
      min-width:30px;
      min-height:30px;
      margin:5px;
      width:85px;
      line-height:30px; }
  .bp3-dark .bp3-file-upload-input{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
    background:rgba(16, 22, 26, 0.3);
    color:#f5f8fa;
    color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input::-webkit-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input::-moz-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input:-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input::-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input::placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input:focus{
      -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-file-upload-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-file-upload-input:disabled, .bp3-dark .bp3-file-upload-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(57, 75, 89, 0.5);
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-file-upload-input::after{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
      background-color:#394b59;
      background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.05)), to(rgba(255, 255, 255, 0)));
      background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.05), rgba(255, 255, 255, 0));
      color:#f5f8fa; }
      .bp3-dark .bp3-file-upload-input::after:hover, .bp3-dark .bp3-file-upload-input::after:active, .bp3-dark .bp3-file-upload-input::after.bp3-active{
        color:#f5f8fa; }
      .bp3-dark .bp3-file-upload-input::after:hover{
        -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
                box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
        background-color:#30404d; }
      .bp3-dark .bp3-file-upload-input::after:active, .bp3-dark .bp3-file-upload-input::after.bp3-active{
        -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
                box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
        background-color:#202b33;
        background-image:none; }
      .bp3-dark .bp3-file-upload-input::after:disabled, .bp3-dark .bp3-file-upload-input::after.bp3-disabled{
        -webkit-box-shadow:none;
                box-shadow:none;
        background-color:rgba(57, 75, 89, 0.5);
        background-image:none;
        color:rgba(167, 182, 194, 0.6); }
        .bp3-dark .bp3-file-upload-input::after:disabled.bp3-active, .bp3-dark .bp3-file-upload-input::after.bp3-disabled.bp3-active{
          background:rgba(57, 75, 89, 0.7); }
      .bp3-dark .bp3-file-upload-input::after .bp3-button-spinner .bp3-spinner-head{
        background:rgba(16, 22, 26, 0.5);
        stroke:#8a9ba8; }
    .bp3-dark .bp3-file-upload-input:hover::after{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
      background-color:#30404d; }
    .bp3-dark .bp3-file-upload-input:active::after{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#202b33;
      background-image:none; }

.bp3-file-upload-input::after{
  -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
          box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1); }
.bp3-form-group{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:vertical;
  -webkit-box-direction:normal;
      -ms-flex-direction:column;
          flex-direction:column;
  margin:0 0 15px; }
  .bp3-form-group label.bp3-label{
    margin-bottom:5px; }
  .bp3-form-group .bp3-control{
    margin-top:7px; }
  .bp3-form-group .bp3-form-helper-text{
    margin-top:5px;
    color:#5c7080;
    font-size:12px; }
  .bp3-form-group.bp3-intent-primary .bp3-form-helper-text{
    color:#106ba3; }
  .bp3-form-group.bp3-intent-success .bp3-form-helper-text{
    color:#0d8050; }
  .bp3-form-group.bp3-intent-warning .bp3-form-helper-text{
    color:#bf7326; }
  .bp3-form-group.bp3-intent-danger .bp3-form-helper-text{
    color:#c23030; }
  .bp3-form-group.bp3-inline{
    -webkit-box-orient:horizontal;
    -webkit-box-direction:normal;
        -ms-flex-direction:row;
            flex-direction:row;
    -webkit-box-align:start;
        -ms-flex-align:start;
            align-items:flex-start; }
    .bp3-form-group.bp3-inline.bp3-large label.bp3-label{
      margin:0 10px 0 0;
      line-height:40px; }
    .bp3-form-group.bp3-inline label.bp3-label{
      margin:0 10px 0 0;
      line-height:30px; }
  .bp3-form-group.bp3-disabled .bp3-label,
  .bp3-form-group.bp3-disabled .bp3-text-muted,
  .bp3-form-group.bp3-disabled .bp3-form-helper-text{
    color:rgba(92, 112, 128, 0.6) !important; }
  .bp3-dark .bp3-form-group.bp3-intent-primary .bp3-form-helper-text{
    color:#48aff0; }
  .bp3-dark .bp3-form-group.bp3-intent-success .bp3-form-helper-text{
    color:#3dcc91; }
  .bp3-dark .bp3-form-group.bp3-intent-warning .bp3-form-helper-text{
    color:#ffb366; }
  .bp3-dark .bp3-form-group.bp3-intent-danger .bp3-form-helper-text{
    color:#ff7373; }
  .bp3-dark .bp3-form-group .bp3-form-helper-text{
    color:#a7b6c2; }
  .bp3-dark .bp3-form-group.bp3-disabled .bp3-label,
  .bp3-dark .bp3-form-group.bp3-disabled .bp3-text-muted,
  .bp3-dark .bp3-form-group.bp3-disabled .bp3-form-helper-text{
    color:rgba(167, 182, 194, 0.6) !important; }
.bp3-input-group{
  display:block;
  position:relative; }
  .bp3-input-group .bp3-input{
    position:relative;
    width:100%; }
    .bp3-input-group .bp3-input:not(:first-child){
      padding-left:30px; }
    .bp3-input-group .bp3-input:not(:last-child){
      padding-right:30px; }
  .bp3-input-group .bp3-input-action,
  .bp3-input-group > .bp3-button,
  .bp3-input-group > .bp3-icon{
    position:absolute;
    top:0; }
    .bp3-input-group .bp3-input-action:first-child,
    .bp3-input-group > .bp3-button:first-child,
    .bp3-input-group > .bp3-icon:first-child{
      left:0; }
    .bp3-input-group .bp3-input-action:last-child,
    .bp3-input-group > .bp3-button:last-child,
    .bp3-input-group > .bp3-icon:last-child{
      right:0; }
  .bp3-input-group .bp3-button{
    min-width:24px;
    min-height:24px;
    margin:3px;
    padding:0 7px; }
    .bp3-input-group .bp3-button:empty{
      padding:0; }
  .bp3-input-group > .bp3-icon{
    z-index:1;
    color:#5c7080; }
    .bp3-input-group > .bp3-icon:empty{
      line-height:1;
      font-family:"Icons16", sans-serif;
      font-size:16px;
      font-weight:400;
      font-style:normal;
      -moz-osx-font-smoothing:grayscale;
      -webkit-font-smoothing:antialiased; }
  .bp3-input-group > .bp3-icon,
  .bp3-input-group .bp3-input-action > .bp3-spinner{
    margin:7px; }
  .bp3-input-group .bp3-tag{
    margin:5px; }
  .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:not(:hover):not(:focus),
  .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:not(:hover):not(:focus){
    color:#5c7080; }
    .bp3-dark .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:not(:hover):not(:focus), .bp3-dark
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:not(:hover):not(:focus){
      color:#a7b6c2; }
    .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:not(:hover):not(:focus) .bp3-icon, .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:not(:hover):not(:focus) .bp3-icon-standard, .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:not(:hover):not(:focus) .bp3-icon-large,
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:not(:hover):not(:focus) .bp3-icon,
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:not(:hover):not(:focus) .bp3-icon-standard,
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:not(:hover):not(:focus) .bp3-icon-large{
      color:#5c7080; }
  .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:disabled,
  .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:disabled{
    color:rgba(92, 112, 128, 0.6) !important; }
    .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:disabled .bp3-icon, .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:disabled .bp3-icon-standard, .bp3-input-group .bp3-input:not(:focus) + .bp3-button.bp3-minimal:disabled .bp3-icon-large,
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:disabled .bp3-icon,
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:disabled .bp3-icon-standard,
    .bp3-input-group .bp3-input:not(:focus) + .bp3-input-action .bp3-button.bp3-minimal:disabled .bp3-icon-large{
      color:rgba(92, 112, 128, 0.6) !important; }
  .bp3-input-group.bp3-disabled{
    cursor:not-allowed; }
    .bp3-input-group.bp3-disabled .bp3-icon{
      color:rgba(92, 112, 128, 0.6); }
  .bp3-input-group.bp3-large .bp3-button{
    min-width:30px;
    min-height:30px;
    margin:5px; }
  .bp3-input-group.bp3-large > .bp3-icon,
  .bp3-input-group.bp3-large .bp3-input-action > .bp3-spinner{
    margin:12px; }
  .bp3-input-group.bp3-large .bp3-input{
    height:40px;
    line-height:40px;
    font-size:16px; }
    .bp3-input-group.bp3-large .bp3-input[type="search"], .bp3-input-group.bp3-large .bp3-input.bp3-round{
      padding:0 15px; }
    .bp3-input-group.bp3-large .bp3-input:not(:first-child){
      padding-left:40px; }
    .bp3-input-group.bp3-large .bp3-input:not(:last-child){
      padding-right:40px; }
  .bp3-input-group.bp3-small .bp3-button{
    min-width:20px;
    min-height:20px;
    margin:2px; }
  .bp3-input-group.bp3-small .bp3-tag{
    min-width:20px;
    min-height:20px;
    margin:2px; }
  .bp3-input-group.bp3-small > .bp3-icon,
  .bp3-input-group.bp3-small .bp3-input-action > .bp3-spinner{
    margin:4px; }
  .bp3-input-group.bp3-small .bp3-input{
    height:24px;
    padding-right:8px;
    padding-left:8px;
    line-height:24px;
    font-size:12px; }
    .bp3-input-group.bp3-small .bp3-input[type="search"], .bp3-input-group.bp3-small .bp3-input.bp3-round{
      padding:0 12px; }
    .bp3-input-group.bp3-small .bp3-input:not(:first-child){
      padding-left:24px; }
    .bp3-input-group.bp3-small .bp3-input:not(:last-child){
      padding-right:24px; }
  .bp3-input-group.bp3-fill{
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto;
    width:100%; }
  .bp3-input-group.bp3-round .bp3-button,
  .bp3-input-group.bp3-round .bp3-input,
  .bp3-input-group.bp3-round .bp3-tag{
    border-radius:30px; }
  .bp3-dark .bp3-input-group .bp3-icon{
    color:#a7b6c2; }
  .bp3-dark .bp3-input-group.bp3-disabled .bp3-icon{
    color:rgba(167, 182, 194, 0.6); }
  .bp3-input-group.bp3-intent-primary .bp3-input{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px #137cbd, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px #137cbd, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-primary .bp3-input:focus{
      -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-primary .bp3-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #137cbd;
              box-shadow:inset 0 0 0 1px #137cbd; }
    .bp3-input-group.bp3-intent-primary .bp3-input:disabled, .bp3-input-group.bp3-intent-primary .bp3-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
  .bp3-input-group.bp3-intent-primary > .bp3-icon{
    color:#106ba3; }
    .bp3-dark .bp3-input-group.bp3-intent-primary > .bp3-icon{
      color:#48aff0; }
  .bp3-input-group.bp3-intent-success .bp3-input{
    -webkit-box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px #0f9960, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px #0f9960, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-success .bp3-input:focus{
      -webkit-box-shadow:0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-success .bp3-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #0f9960;
              box-shadow:inset 0 0 0 1px #0f9960; }
    .bp3-input-group.bp3-intent-success .bp3-input:disabled, .bp3-input-group.bp3-intent-success .bp3-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
  .bp3-input-group.bp3-intent-success > .bp3-icon{
    color:#0d8050; }
    .bp3-dark .bp3-input-group.bp3-intent-success > .bp3-icon{
      color:#3dcc91; }
  .bp3-input-group.bp3-intent-warning .bp3-input{
    -webkit-box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px #d9822b, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px #d9822b, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-warning .bp3-input:focus{
      -webkit-box-shadow:0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-warning .bp3-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #d9822b;
              box-shadow:inset 0 0 0 1px #d9822b; }
    .bp3-input-group.bp3-intent-warning .bp3-input:disabled, .bp3-input-group.bp3-intent-warning .bp3-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
  .bp3-input-group.bp3-intent-warning > .bp3-icon{
    color:#bf7326; }
    .bp3-dark .bp3-input-group.bp3-intent-warning > .bp3-icon{
      color:#ffb366; }
  .bp3-input-group.bp3-intent-danger .bp3-input{
    -webkit-box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px #db3737, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px #db3737, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-danger .bp3-input:focus{
      -webkit-box-shadow:0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input-group.bp3-intent-danger .bp3-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #db3737;
              box-shadow:inset 0 0 0 1px #db3737; }
    .bp3-input-group.bp3-intent-danger .bp3-input:disabled, .bp3-input-group.bp3-intent-danger .bp3-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
  .bp3-input-group.bp3-intent-danger > .bp3-icon{
    color:#c23030; }
    .bp3-dark .bp3-input-group.bp3-intent-danger > .bp3-icon{
      color:#ff7373; }
.bp3-input{
  outline:none;
  border:none;
  border-radius:3px;
  -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
  background:#ffffff;
  height:30px;
  padding:0 10px;
  vertical-align:middle;
  line-height:30px;
  color:#182026;
  font-size:14px;
  font-weight:400;
  -webkit-transition:-webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:-webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-box-shadow 100ms cubic-bezier(0.4, 1, 0.75, 0.9);
  -webkit-appearance:none;
     -moz-appearance:none;
          appearance:none; }
  .bp3-input::-webkit-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input::-moz-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input:-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input::-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input::placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input:focus, .bp3-input.bp3-active{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-input[type="search"], .bp3-input.bp3-round{
    border-radius:30px;
    -webkit-box-sizing:border-box;
            box-sizing:border-box;
    padding-left:10px; }
  .bp3-input[readonly]{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.15);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.15); }
  .bp3-input:disabled, .bp3-input.bp3-disabled{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(206, 217, 224, 0.5);
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6);
    resize:none; }
  .bp3-input.bp3-large{
    height:40px;
    line-height:40px;
    font-size:16px; }
    .bp3-input.bp3-large[type="search"], .bp3-input.bp3-large.bp3-round{
      padding:0 15px; }
  .bp3-input.bp3-small{
    height:24px;
    padding-right:8px;
    padding-left:8px;
    line-height:24px;
    font-size:12px; }
    .bp3-input.bp3-small[type="search"], .bp3-input.bp3-small.bp3-round{
      padding:0 12px; }
  .bp3-input.bp3-fill{
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto;
    width:100%; }
  .bp3-dark .bp3-input{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
    background:rgba(16, 22, 26, 0.3);
    color:#f5f8fa; }
    .bp3-dark .bp3-input::-webkit-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-input::-moz-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-input:-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-input::-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-input::placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-input:focus{
      -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-input:disabled, .bp3-dark .bp3-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(57, 75, 89, 0.5);
      color:rgba(167, 182, 194, 0.6); }
  .bp3-input.bp3-intent-primary{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px #137cbd, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px #137cbd, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-primary:focus{
      -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-primary[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #137cbd;
              box-shadow:inset 0 0 0 1px #137cbd; }
    .bp3-input.bp3-intent-primary:disabled, .bp3-input.bp3-intent-primary.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
    .bp3-dark .bp3-input.bp3-intent-primary{
      -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px #137cbd, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px #137cbd, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-primary:focus{
        -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
                box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-primary[readonly]{
        -webkit-box-shadow:inset 0 0 0 1px #137cbd;
                box-shadow:inset 0 0 0 1px #137cbd; }
      .bp3-dark .bp3-input.bp3-intent-primary:disabled, .bp3-dark .bp3-input.bp3-intent-primary.bp3-disabled{
        -webkit-box-shadow:none;
                box-shadow:none; }
  .bp3-input.bp3-intent-success{
    -webkit-box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px #0f9960, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px #0f9960, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-success:focus{
      -webkit-box-shadow:0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-success[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #0f9960;
              box-shadow:inset 0 0 0 1px #0f9960; }
    .bp3-input.bp3-intent-success:disabled, .bp3-input.bp3-intent-success.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
    .bp3-dark .bp3-input.bp3-intent-success{
      -webkit-box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px #0f9960, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), 0 0 0 0 rgba(15, 153, 96, 0), inset 0 0 0 1px #0f9960, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-success:focus{
        -webkit-box-shadow:0 0 0 1px #0f9960, 0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
                box-shadow:0 0 0 1px #0f9960, 0 0 0 1px #0f9960, 0 0 0 3px rgba(15, 153, 96, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-success[readonly]{
        -webkit-box-shadow:inset 0 0 0 1px #0f9960;
                box-shadow:inset 0 0 0 1px #0f9960; }
      .bp3-dark .bp3-input.bp3-intent-success:disabled, .bp3-dark .bp3-input.bp3-intent-success.bp3-disabled{
        -webkit-box-shadow:none;
                box-shadow:none; }
  .bp3-input.bp3-intent-warning{
    -webkit-box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px #d9822b, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px #d9822b, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-warning:focus{
      -webkit-box-shadow:0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-warning[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #d9822b;
              box-shadow:inset 0 0 0 1px #d9822b; }
    .bp3-input.bp3-intent-warning:disabled, .bp3-input.bp3-intent-warning.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
    .bp3-dark .bp3-input.bp3-intent-warning{
      -webkit-box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px #d9822b, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), 0 0 0 0 rgba(217, 130, 43, 0), inset 0 0 0 1px #d9822b, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-warning:focus{
        -webkit-box-shadow:0 0 0 1px #d9822b, 0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
                box-shadow:0 0 0 1px #d9822b, 0 0 0 1px #d9822b, 0 0 0 3px rgba(217, 130, 43, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-warning[readonly]{
        -webkit-box-shadow:inset 0 0 0 1px #d9822b;
                box-shadow:inset 0 0 0 1px #d9822b; }
      .bp3-dark .bp3-input.bp3-intent-warning:disabled, .bp3-dark .bp3-input.bp3-intent-warning.bp3-disabled{
        -webkit-box-shadow:none;
                box-shadow:none; }
  .bp3-input.bp3-intent-danger{
    -webkit-box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px #db3737, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px #db3737, inset 0 0 0 1px rgba(16, 22, 26, 0.15), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-danger:focus{
      -webkit-box-shadow:0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-input.bp3-intent-danger[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px #db3737;
              box-shadow:inset 0 0 0 1px #db3737; }
    .bp3-input.bp3-intent-danger:disabled, .bp3-input.bp3-intent-danger.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none; }
    .bp3-dark .bp3-input.bp3-intent-danger{
      -webkit-box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px #db3737, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), 0 0 0 0 rgba(219, 55, 55, 0), inset 0 0 0 1px #db3737, inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-danger:focus{
        -webkit-box-shadow:0 0 0 1px #db3737, 0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
                box-shadow:0 0 0 1px #db3737, 0 0 0 1px #db3737, 0 0 0 3px rgba(219, 55, 55, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
      .bp3-dark .bp3-input.bp3-intent-danger[readonly]{
        -webkit-box-shadow:inset 0 0 0 1px #db3737;
                box-shadow:inset 0 0 0 1px #db3737; }
      .bp3-dark .bp3-input.bp3-intent-danger:disabled, .bp3-dark .bp3-input.bp3-intent-danger.bp3-disabled{
        -webkit-box-shadow:none;
                box-shadow:none; }
  .bp3-input::-ms-clear{
    display:none; }
textarea.bp3-input{
  max-width:100%;
  padding:10px; }
  textarea.bp3-input, textarea.bp3-input.bp3-large, textarea.bp3-input.bp3-small{
    height:auto;
    line-height:inherit; }
  textarea.bp3-input.bp3-small{
    padding:8px; }
  .bp3-dark textarea.bp3-input{
    -webkit-box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), 0 0 0 0 rgba(19, 124, 189, 0), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
    background:rgba(16, 22, 26, 0.3);
    color:#f5f8fa; }
    .bp3-dark textarea.bp3-input::-webkit-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark textarea.bp3-input::-moz-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark textarea.bp3-input:-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark textarea.bp3-input::-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark textarea.bp3-input::placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark textarea.bp3-input:focus{
      -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark textarea.bp3-input[readonly]{
      -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark textarea.bp3-input:disabled, .bp3-dark textarea.bp3-input.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:rgba(57, 75, 89, 0.5);
      color:rgba(167, 182, 194, 0.6); }
label.bp3-label{
  display:block;
  margin-top:0;
  margin-bottom:15px; }
  label.bp3-label .bp3-html-select,
  label.bp3-label .bp3-input,
  label.bp3-label .bp3-select,
  label.bp3-label .bp3-slider,
  label.bp3-label .bp3-popover-wrapper{
    display:block;
    margin-top:5px;
    text-transform:none; }
  label.bp3-label .bp3-button-group{
    margin-top:5px; }
  label.bp3-label .bp3-select select,
  label.bp3-label .bp3-html-select select{
    width:100%;
    vertical-align:top;
    font-weight:400; }
  label.bp3-label.bp3-disabled,
  label.bp3-label.bp3-disabled .bp3-text-muted{
    color:rgba(92, 112, 128, 0.6); }
  label.bp3-label.bp3-inline{
    line-height:30px; }
    label.bp3-label.bp3-inline .bp3-html-select,
    label.bp3-label.bp3-inline .bp3-input,
    label.bp3-label.bp3-inline .bp3-input-group,
    label.bp3-label.bp3-inline .bp3-select,
    label.bp3-label.bp3-inline .bp3-popover-wrapper{
      display:inline-block;
      margin:0 0 0 5px;
      vertical-align:top; }
    label.bp3-label.bp3-inline .bp3-button-group{
      margin:0 0 0 5px; }
    label.bp3-label.bp3-inline .bp3-input-group .bp3-input{
      margin-left:0; }
    label.bp3-label.bp3-inline.bp3-large{
      line-height:40px; }
  label.bp3-label:not(.bp3-inline) .bp3-popover-target{
    display:block; }
  .bp3-dark label.bp3-label{
    color:#f5f8fa; }
    .bp3-dark label.bp3-label.bp3-disabled,
    .bp3-dark label.bp3-label.bp3-disabled .bp3-text-muted{
      color:rgba(167, 182, 194, 0.6); }
.bp3-numeric-input .bp3-button-group.bp3-vertical > .bp3-button{
  -webkit-box-flex:1;
      -ms-flex:1 1 14px;
          flex:1 1 14px;
  width:30px;
  min-height:0;
  padding:0; }
  .bp3-numeric-input .bp3-button-group.bp3-vertical > .bp3-button:first-child{
    border-radius:0 3px 0 0; }
  .bp3-numeric-input .bp3-button-group.bp3-vertical > .bp3-button:last-child{
    border-radius:0 0 3px 0; }

.bp3-numeric-input .bp3-button-group.bp3-vertical:first-child > .bp3-button:first-child{
  border-radius:3px 0 0 0; }

.bp3-numeric-input .bp3-button-group.bp3-vertical:first-child > .bp3-button:last-child{
  border-radius:0 0 0 3px; }

.bp3-numeric-input.bp3-large .bp3-button-group.bp3-vertical > .bp3-button{
  width:40px; }

form{
  display:block; }
.bp3-html-select select,
.bp3-select select{
  display:-webkit-inline-box;
  display:-ms-inline-flexbox;
  display:inline-flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:center;
      -ms-flex-pack:center;
          justify-content:center;
  border:none;
  border-radius:3px;
  cursor:pointer;
  padding:5px 10px;
  vertical-align:middle;
  text-align:left;
  font-size:14px;
  -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
          box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
  background-color:#f5f8fa;
  background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.8)), to(rgba(255, 255, 255, 0)));
  background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.8), rgba(255, 255, 255, 0));
  color:#182026;
  border-radius:3px;
  width:100%;
  height:30px;
  padding:0 25px 0 10px;
  -moz-appearance:none;
  -webkit-appearance:none; }
  .bp3-html-select select > *, .bp3-select select > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-html-select select > .bp3-fill, .bp3-select select > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-html-select select::before,
  .bp3-select select::before, .bp3-html-select select > *, .bp3-select select > *{
    margin-right:7px; }
  .bp3-html-select select:empty::before,
  .bp3-select select:empty::before,
  .bp3-html-select select > :last-child,
  .bp3-select select > :last-child{
    margin-right:0; }
  .bp3-html-select select:hover,
  .bp3-select select:hover{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-clip:padding-box;
    background-color:#ebf1f5; }
  .bp3-html-select select:active,
  .bp3-select select:active, .bp3-html-select select.bp3-active,
  .bp3-select select.bp3-active{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#d8e1e8;
    background-image:none; }
  .bp3-html-select select:disabled,
  .bp3-select select:disabled, .bp3-html-select select.bp3-disabled,
  .bp3-select select.bp3-disabled{
    outline:none;
    -webkit-box-shadow:none;
            box-shadow:none;
    background-color:rgba(206, 217, 224, 0.5);
    background-image:none;
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6); }
    .bp3-html-select select:disabled.bp3-active,
    .bp3-select select:disabled.bp3-active, .bp3-html-select select:disabled.bp3-active:hover,
    .bp3-select select:disabled.bp3-active:hover, .bp3-html-select select.bp3-disabled.bp3-active,
    .bp3-select select.bp3-disabled.bp3-active, .bp3-html-select select.bp3-disabled.bp3-active:hover,
    .bp3-select select.bp3-disabled.bp3-active:hover{
      background:rgba(206, 217, 224, 0.7); }

.bp3-html-select.bp3-minimal select,
.bp3-select.bp3-minimal select{
  -webkit-box-shadow:none;
          box-shadow:none;
  background:none; }
  .bp3-html-select.bp3-minimal select:hover,
  .bp3-select.bp3-minimal select:hover{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(167, 182, 194, 0.3);
    text-decoration:none;
    color:#182026; }
  .bp3-html-select.bp3-minimal select:active,
  .bp3-select.bp3-minimal select:active, .bp3-html-select.bp3-minimal select.bp3-active,
  .bp3-select.bp3-minimal select.bp3-active{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:rgba(115, 134, 148, 0.3);
    color:#182026; }
  .bp3-html-select.bp3-minimal select:disabled,
  .bp3-select.bp3-minimal select:disabled, .bp3-html-select.bp3-minimal select:disabled:hover,
  .bp3-select.bp3-minimal select:disabled:hover, .bp3-html-select.bp3-minimal select.bp3-disabled,
  .bp3-select.bp3-minimal select.bp3-disabled, .bp3-html-select.bp3-minimal select.bp3-disabled:hover,
  .bp3-select.bp3-minimal select.bp3-disabled:hover{
    background:none;
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6); }
    .bp3-html-select.bp3-minimal select:disabled.bp3-active,
    .bp3-select.bp3-minimal select:disabled.bp3-active, .bp3-html-select.bp3-minimal select:disabled:hover.bp3-active,
    .bp3-select.bp3-minimal select:disabled:hover.bp3-active, .bp3-html-select.bp3-minimal select.bp3-disabled.bp3-active,
    .bp3-select.bp3-minimal select.bp3-disabled.bp3-active, .bp3-html-select.bp3-minimal select.bp3-disabled:hover.bp3-active,
    .bp3-select.bp3-minimal select.bp3-disabled:hover.bp3-active{
      background:rgba(115, 134, 148, 0.3); }
  .bp3-dark .bp3-html-select.bp3-minimal select, .bp3-html-select.bp3-minimal .bp3-dark select,
  .bp3-dark .bp3-select.bp3-minimal select, .bp3-select.bp3-minimal .bp3-dark select{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:none;
    color:inherit; }
    .bp3-dark .bp3-html-select.bp3-minimal select:hover, .bp3-html-select.bp3-minimal .bp3-dark select:hover,
    .bp3-dark .bp3-select.bp3-minimal select:hover, .bp3-select.bp3-minimal .bp3-dark select:hover, .bp3-dark .bp3-html-select.bp3-minimal select:active, .bp3-html-select.bp3-minimal .bp3-dark select:active,
    .bp3-dark .bp3-select.bp3-minimal select:active, .bp3-select.bp3-minimal .bp3-dark select:active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-active,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none; }
    .bp3-dark .bp3-html-select.bp3-minimal select:hover, .bp3-html-select.bp3-minimal .bp3-dark select:hover,
    .bp3-dark .bp3-select.bp3-minimal select:hover, .bp3-select.bp3-minimal .bp3-dark select:hover{
      background:rgba(138, 155, 168, 0.15); }
    .bp3-dark .bp3-html-select.bp3-minimal select:active, .bp3-html-select.bp3-minimal .bp3-dark select:active,
    .bp3-dark .bp3-select.bp3-minimal select:active, .bp3-select.bp3-minimal .bp3-dark select:active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-active,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-active{
      background:rgba(138, 155, 168, 0.3);
      color:#f5f8fa; }
    .bp3-dark .bp3-html-select.bp3-minimal select:disabled, .bp3-html-select.bp3-minimal .bp3-dark select:disabled,
    .bp3-dark .bp3-select.bp3-minimal select:disabled, .bp3-select.bp3-minimal .bp3-dark select:disabled, .bp3-dark .bp3-html-select.bp3-minimal select:disabled:hover, .bp3-html-select.bp3-minimal .bp3-dark select:disabled:hover,
    .bp3-dark .bp3-select.bp3-minimal select:disabled:hover, .bp3-select.bp3-minimal .bp3-dark select:disabled:hover, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-disabled,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-disabled, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-disabled:hover, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-disabled:hover,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-disabled:hover, .bp3-select.bp3-minimal .bp3-dark select.bp3-disabled:hover{
      background:none;
      cursor:not-allowed;
      color:rgba(167, 182, 194, 0.6); }
      .bp3-dark .bp3-html-select.bp3-minimal select:disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select:disabled.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select:disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select:disabled.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select:disabled:hover.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select:disabled:hover.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select:disabled:hover.bp3-active, .bp3-select.bp3-minimal .bp3-dark select:disabled:hover.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-disabled.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-disabled.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-disabled:hover.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-disabled:hover.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-disabled:hover.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-disabled:hover.bp3-active{
        background:rgba(138, 155, 168, 0.3); }
  .bp3-html-select.bp3-minimal select.bp3-intent-primary,
  .bp3-select.bp3-minimal select.bp3-intent-primary{
    color:#106ba3; }
    .bp3-html-select.bp3-minimal select.bp3-intent-primary:hover,
    .bp3-select.bp3-minimal select.bp3-intent-primary:hover, .bp3-html-select.bp3-minimal select.bp3-intent-primary:active,
    .bp3-select.bp3-minimal select.bp3-intent-primary:active, .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none;
      color:#106ba3; }
    .bp3-html-select.bp3-minimal select.bp3-intent-primary:hover,
    .bp3-select.bp3-minimal select.bp3-intent-primary:hover{
      background:rgba(19, 124, 189, 0.15);
      color:#106ba3; }
    .bp3-html-select.bp3-minimal select.bp3-intent-primary:active,
    .bp3-select.bp3-minimal select.bp3-intent-primary:active, .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-active{
      background:rgba(19, 124, 189, 0.3);
      color:#106ba3; }
    .bp3-html-select.bp3-minimal select.bp3-intent-primary:disabled,
    .bp3-select.bp3-minimal select.bp3-intent-primary:disabled, .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-disabled,
    .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-disabled{
      background:none;
      color:rgba(16, 107, 163, 0.5); }
      .bp3-html-select.bp3-minimal select.bp3-intent-primary:disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-primary:disabled.bp3-active, .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-disabled.bp3-active{
        background:rgba(19, 124, 189, 0.3); }
    .bp3-html-select.bp3-minimal select.bp3-intent-primary .bp3-button-spinner .bp3-spinner-head, .bp3-select.bp3-minimal select.bp3-intent-primary .bp3-button-spinner .bp3-spinner-head{
      stroke:#106ba3; }
    .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary{
      color:#48aff0; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary:hover, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary:hover,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary:hover, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary:hover{
        background:rgba(19, 124, 189, 0.2);
        color:#48aff0; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary:active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary:active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary:active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary:active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary.bp3-active{
        background:rgba(19, 124, 189, 0.3);
        color:#48aff0; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary:disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary:disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary:disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary:disabled, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary.bp3-disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary.bp3-disabled{
        background:none;
        color:rgba(72, 175, 240, 0.5); }
        .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary:disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary:disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary:disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary:disabled.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-primary.bp3-disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-primary.bp3-disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-primary.bp3-disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-primary.bp3-disabled.bp3-active{
          background:rgba(19, 124, 189, 0.3); }
  .bp3-html-select.bp3-minimal select.bp3-intent-success,
  .bp3-select.bp3-minimal select.bp3-intent-success{
    color:#0d8050; }
    .bp3-html-select.bp3-minimal select.bp3-intent-success:hover,
    .bp3-select.bp3-minimal select.bp3-intent-success:hover, .bp3-html-select.bp3-minimal select.bp3-intent-success:active,
    .bp3-select.bp3-minimal select.bp3-intent-success:active, .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-success.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none;
      color:#0d8050; }
    .bp3-html-select.bp3-minimal select.bp3-intent-success:hover,
    .bp3-select.bp3-minimal select.bp3-intent-success:hover{
      background:rgba(15, 153, 96, 0.15);
      color:#0d8050; }
    .bp3-html-select.bp3-minimal select.bp3-intent-success:active,
    .bp3-select.bp3-minimal select.bp3-intent-success:active, .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-success.bp3-active{
      background:rgba(15, 153, 96, 0.3);
      color:#0d8050; }
    .bp3-html-select.bp3-minimal select.bp3-intent-success:disabled,
    .bp3-select.bp3-minimal select.bp3-intent-success:disabled, .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-disabled,
    .bp3-select.bp3-minimal select.bp3-intent-success.bp3-disabled{
      background:none;
      color:rgba(13, 128, 80, 0.5); }
      .bp3-html-select.bp3-minimal select.bp3-intent-success:disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-success:disabled.bp3-active, .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-success.bp3-disabled.bp3-active{
        background:rgba(15, 153, 96, 0.3); }
    .bp3-html-select.bp3-minimal select.bp3-intent-success .bp3-button-spinner .bp3-spinner-head, .bp3-select.bp3-minimal select.bp3-intent-success .bp3-button-spinner .bp3-spinner-head{
      stroke:#0d8050; }
    .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success{
      color:#3dcc91; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success:hover, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success:hover,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success:hover, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success:hover{
        background:rgba(15, 153, 96, 0.2);
        color:#3dcc91; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success:active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success:active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success:active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success:active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success.bp3-active{
        background:rgba(15, 153, 96, 0.3);
        color:#3dcc91; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success:disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success:disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success:disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success:disabled, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success.bp3-disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success.bp3-disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success.bp3-disabled{
        background:none;
        color:rgba(61, 204, 145, 0.5); }
        .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success:disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success:disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success:disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success:disabled.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-success.bp3-disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-success.bp3-disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-success.bp3-disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-success.bp3-disabled.bp3-active{
          background:rgba(15, 153, 96, 0.3); }
  .bp3-html-select.bp3-minimal select.bp3-intent-warning,
  .bp3-select.bp3-minimal select.bp3-intent-warning{
    color:#bf7326; }
    .bp3-html-select.bp3-minimal select.bp3-intent-warning:hover,
    .bp3-select.bp3-minimal select.bp3-intent-warning:hover, .bp3-html-select.bp3-minimal select.bp3-intent-warning:active,
    .bp3-select.bp3-minimal select.bp3-intent-warning:active, .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none;
      color:#bf7326; }
    .bp3-html-select.bp3-minimal select.bp3-intent-warning:hover,
    .bp3-select.bp3-minimal select.bp3-intent-warning:hover{
      background:rgba(217, 130, 43, 0.15);
      color:#bf7326; }
    .bp3-html-select.bp3-minimal select.bp3-intent-warning:active,
    .bp3-select.bp3-minimal select.bp3-intent-warning:active, .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-active{
      background:rgba(217, 130, 43, 0.3);
      color:#bf7326; }
    .bp3-html-select.bp3-minimal select.bp3-intent-warning:disabled,
    .bp3-select.bp3-minimal select.bp3-intent-warning:disabled, .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-disabled,
    .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-disabled{
      background:none;
      color:rgba(191, 115, 38, 0.5); }
      .bp3-html-select.bp3-minimal select.bp3-intent-warning:disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-warning:disabled.bp3-active, .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-disabled.bp3-active{
        background:rgba(217, 130, 43, 0.3); }
    .bp3-html-select.bp3-minimal select.bp3-intent-warning .bp3-button-spinner .bp3-spinner-head, .bp3-select.bp3-minimal select.bp3-intent-warning .bp3-button-spinner .bp3-spinner-head{
      stroke:#bf7326; }
    .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning{
      color:#ffb366; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning:hover, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning:hover,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning:hover, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning:hover{
        background:rgba(217, 130, 43, 0.2);
        color:#ffb366; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning:active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning:active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning:active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning:active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning.bp3-active{
        background:rgba(217, 130, 43, 0.3);
        color:#ffb366; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning:disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning:disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning:disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning:disabled, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning.bp3-disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning.bp3-disabled{
        background:none;
        color:rgba(255, 179, 102, 0.5); }
        .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning:disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning:disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning:disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning:disabled.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-warning.bp3-disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-warning.bp3-disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-warning.bp3-disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-warning.bp3-disabled.bp3-active{
          background:rgba(217, 130, 43, 0.3); }
  .bp3-html-select.bp3-minimal select.bp3-intent-danger,
  .bp3-select.bp3-minimal select.bp3-intent-danger{
    color:#c23030; }
    .bp3-html-select.bp3-minimal select.bp3-intent-danger:hover,
    .bp3-select.bp3-minimal select.bp3-intent-danger:hover, .bp3-html-select.bp3-minimal select.bp3-intent-danger:active,
    .bp3-select.bp3-minimal select.bp3-intent-danger:active, .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-active{
      -webkit-box-shadow:none;
              box-shadow:none;
      background:none;
      color:#c23030; }
    .bp3-html-select.bp3-minimal select.bp3-intent-danger:hover,
    .bp3-select.bp3-minimal select.bp3-intent-danger:hover{
      background:rgba(219, 55, 55, 0.15);
      color:#c23030; }
    .bp3-html-select.bp3-minimal select.bp3-intent-danger:active,
    .bp3-select.bp3-minimal select.bp3-intent-danger:active, .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-active,
    .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-active{
      background:rgba(219, 55, 55, 0.3);
      color:#c23030; }
    .bp3-html-select.bp3-minimal select.bp3-intent-danger:disabled,
    .bp3-select.bp3-minimal select.bp3-intent-danger:disabled, .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-disabled,
    .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-disabled{
      background:none;
      color:rgba(194, 48, 48, 0.5); }
      .bp3-html-select.bp3-minimal select.bp3-intent-danger:disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-danger:disabled.bp3-active, .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-disabled.bp3-active,
      .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-disabled.bp3-active{
        background:rgba(219, 55, 55, 0.3); }
    .bp3-html-select.bp3-minimal select.bp3-intent-danger .bp3-button-spinner .bp3-spinner-head, .bp3-select.bp3-minimal select.bp3-intent-danger .bp3-button-spinner .bp3-spinner-head{
      stroke:#c23030; }
    .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger,
    .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger{
      color:#ff7373; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger:hover, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger:hover,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger:hover, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger:hover{
        background:rgba(219, 55, 55, 0.2);
        color:#ff7373; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger:active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger:active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger:active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger:active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger.bp3-active,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger.bp3-active{
        background:rgba(219, 55, 55, 0.3);
        color:#ff7373; }
      .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger:disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger:disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger:disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger:disabled, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-disabled, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger.bp3-disabled,
      .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-disabled, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger.bp3-disabled{
        background:none;
        color:rgba(255, 115, 115, 0.5); }
        .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger:disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger:disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger:disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger:disabled.bp3-active, .bp3-dark .bp3-html-select.bp3-minimal select.bp3-intent-danger.bp3-disabled.bp3-active, .bp3-html-select.bp3-minimal .bp3-dark select.bp3-intent-danger.bp3-disabled.bp3-active,
        .bp3-dark .bp3-select.bp3-minimal select.bp3-intent-danger.bp3-disabled.bp3-active, .bp3-select.bp3-minimal .bp3-dark select.bp3-intent-danger.bp3-disabled.bp3-active{
          background:rgba(219, 55, 55, 0.3); }

.bp3-html-select.bp3-large select,
.bp3-select.bp3-large select{
  height:40px;
  padding-right:35px;
  font-size:16px; }

.bp3-dark .bp3-html-select select, .bp3-dark .bp3-select select{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
  background-color:#394b59;
  background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.05)), to(rgba(255, 255, 255, 0)));
  background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.05), rgba(255, 255, 255, 0));
  color:#f5f8fa; }
  .bp3-dark .bp3-html-select select:hover, .bp3-dark .bp3-select select:hover, .bp3-dark .bp3-html-select select:active, .bp3-dark .bp3-select select:active, .bp3-dark .bp3-html-select select.bp3-active, .bp3-dark .bp3-select select.bp3-active{
    color:#f5f8fa; }
  .bp3-dark .bp3-html-select select:hover, .bp3-dark .bp3-select select:hover{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
    background-color:#30404d; }
  .bp3-dark .bp3-html-select select:active, .bp3-dark .bp3-select select:active, .bp3-dark .bp3-html-select select.bp3-active, .bp3-dark .bp3-select select.bp3-active{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#202b33;
    background-image:none; }
  .bp3-dark .bp3-html-select select:disabled, .bp3-dark .bp3-select select:disabled, .bp3-dark .bp3-html-select select.bp3-disabled, .bp3-dark .bp3-select select.bp3-disabled{
    -webkit-box-shadow:none;
            box-shadow:none;
    background-color:rgba(57, 75, 89, 0.5);
    background-image:none;
    color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-html-select select:disabled.bp3-active, .bp3-dark .bp3-select select:disabled.bp3-active, .bp3-dark .bp3-html-select select.bp3-disabled.bp3-active, .bp3-dark .bp3-select select.bp3-disabled.bp3-active{
      background:rgba(57, 75, 89, 0.7); }
  .bp3-dark .bp3-html-select select .bp3-button-spinner .bp3-spinner-head, .bp3-dark .bp3-select select .bp3-button-spinner .bp3-spinner-head{
    background:rgba(16, 22, 26, 0.5);
    stroke:#8a9ba8; }

.bp3-html-select select:disabled,
.bp3-select select:disabled{
  -webkit-box-shadow:none;
          box-shadow:none;
  background-color:rgba(206, 217, 224, 0.5);
  cursor:not-allowed;
  color:rgba(92, 112, 128, 0.6); }

.bp3-html-select .bp3-icon,
.bp3-select .bp3-icon, .bp3-select::after{
  position:absolute;
  top:7px;
  right:7px;
  color:#5c7080;
  pointer-events:none; }
  .bp3-html-select .bp3-disabled.bp3-icon,
  .bp3-select .bp3-disabled.bp3-icon, .bp3-disabled.bp3-select::after{
    color:rgba(92, 112, 128, 0.6); }
.bp3-html-select,
.bp3-select{
  display:inline-block;
  position:relative;
  vertical-align:middle;
  letter-spacing:normal; }
  .bp3-html-select select::-ms-expand,
  .bp3-select select::-ms-expand{
    display:none; }
  .bp3-html-select .bp3-icon,
  .bp3-select .bp3-icon{
    color:#5c7080; }
    .bp3-html-select .bp3-icon:hover,
    .bp3-select .bp3-icon:hover{
      color:#182026; }
    .bp3-dark .bp3-html-select .bp3-icon, .bp3-dark
    .bp3-select .bp3-icon{
      color:#a7b6c2; }
      .bp3-dark .bp3-html-select .bp3-icon:hover, .bp3-dark
      .bp3-select .bp3-icon:hover{
        color:#f5f8fa; }
  .bp3-html-select.bp3-large::after,
  .bp3-html-select.bp3-large .bp3-icon,
  .bp3-select.bp3-large::after,
  .bp3-select.bp3-large .bp3-icon{
    top:12px;
    right:12px; }
  .bp3-html-select.bp3-fill,
  .bp3-html-select.bp3-fill select,
  .bp3-select.bp3-fill,
  .bp3-select.bp3-fill select{
    width:100%; }
  .bp3-dark .bp3-html-select option, .bp3-dark
  .bp3-select option{
    background-color:#30404d;
    color:#f5f8fa; }
  .bp3-dark .bp3-html-select::after, .bp3-dark
  .bp3-select::after{
    color:#a7b6c2; }

.bp3-select::after{
  line-height:1;
  font-family:"Icons16", sans-serif;
  font-size:16px;
  font-weight:400;
  font-style:normal;
  -moz-osx-font-smoothing:grayscale;
  -webkit-font-smoothing:antialiased;
  content:""; }
.bp3-running-text table, table.bp3-html-table{
  border-spacing:0;
  font-size:14px; }
  .bp3-running-text table th, table.bp3-html-table th,
  .bp3-running-text table td,
  table.bp3-html-table td{
    padding:11px;
    vertical-align:top;
    text-align:left; }
  .bp3-running-text table th, table.bp3-html-table th{
    color:#182026;
    font-weight:600; }
  
  .bp3-running-text table td,
  table.bp3-html-table td{
    color:#182026; }
  .bp3-running-text table tbody tr:first-child th, table.bp3-html-table tbody tr:first-child th,
  .bp3-running-text table tbody tr:first-child td,
  table.bp3-html-table tbody tr:first-child td{
    -webkit-box-shadow:inset 0 1px 0 0 rgba(16, 22, 26, 0.15);
            box-shadow:inset 0 1px 0 0 rgba(16, 22, 26, 0.15); }
  .bp3-dark .bp3-running-text table th, .bp3-running-text .bp3-dark table th, .bp3-dark table.bp3-html-table th{
    color:#f5f8fa; }
  .bp3-dark .bp3-running-text table td, .bp3-running-text .bp3-dark table td, .bp3-dark table.bp3-html-table td{
    color:#f5f8fa; }
  .bp3-dark .bp3-running-text table tbody tr:first-child th, .bp3-running-text .bp3-dark table tbody tr:first-child th, .bp3-dark table.bp3-html-table tbody tr:first-child th,
  .bp3-dark .bp3-running-text table tbody tr:first-child td,
  .bp3-running-text .bp3-dark table tbody tr:first-child td,
  .bp3-dark table.bp3-html-table tbody tr:first-child td{
    -webkit-box-shadow:inset 0 1px 0 0 rgba(255, 255, 255, 0.15);
            box-shadow:inset 0 1px 0 0 rgba(255, 255, 255, 0.15); }

table.bp3-html-table.bp3-html-table-condensed th,
table.bp3-html-table.bp3-html-table-condensed td, table.bp3-html-table.bp3-small th,
table.bp3-html-table.bp3-small td{
  padding-top:6px;
  padding-bottom:6px; }

table.bp3-html-table.bp3-html-table-striped tbody tr:nth-child(odd) td{
  background:rgba(191, 204, 214, 0.15); }

table.bp3-html-table.bp3-html-table-bordered th:not(:first-child){
  -webkit-box-shadow:inset 1px 0 0 0 rgba(16, 22, 26, 0.15);
          box-shadow:inset 1px 0 0 0 rgba(16, 22, 26, 0.15); }

table.bp3-html-table.bp3-html-table-bordered tbody tr td{
  -webkit-box-shadow:inset 0 1px 0 0 rgba(16, 22, 26, 0.15);
          box-shadow:inset 0 1px 0 0 rgba(16, 22, 26, 0.15); }
  table.bp3-html-table.bp3-html-table-bordered tbody tr td:not(:first-child){
    -webkit-box-shadow:inset 1px 1px 0 0 rgba(16, 22, 26, 0.15);
            box-shadow:inset 1px 1px 0 0 rgba(16, 22, 26, 0.15); }

table.bp3-html-table.bp3-html-table-bordered.bp3-html-table-striped tbody tr:not(:first-child) td{
  -webkit-box-shadow:none;
          box-shadow:none; }
  table.bp3-html-table.bp3-html-table-bordered.bp3-html-table-striped tbody tr:not(:first-child) td:not(:first-child){
    -webkit-box-shadow:inset 1px 0 0 0 rgba(16, 22, 26, 0.15);
            box-shadow:inset 1px 0 0 0 rgba(16, 22, 26, 0.15); }

table.bp3-html-table.bp3-interactive tbody tr:hover td{
  background-color:rgba(191, 204, 214, 0.3);
  cursor:pointer; }

table.bp3-html-table.bp3-interactive tbody tr:active td{
  background-color:rgba(191, 204, 214, 0.4); }

.bp3-dark table.bp3-html-table.bp3-html-table-striped tbody tr:nth-child(odd) td{
  background:rgba(92, 112, 128, 0.15); }

.bp3-dark table.bp3-html-table.bp3-html-table-bordered th:not(:first-child){
  -webkit-box-shadow:inset 1px 0 0 0 rgba(255, 255, 255, 0.15);
          box-shadow:inset 1px 0 0 0 rgba(255, 255, 255, 0.15); }

.bp3-dark table.bp3-html-table.bp3-html-table-bordered tbody tr td{
  -webkit-box-shadow:inset 0 1px 0 0 rgba(255, 255, 255, 0.15);
          box-shadow:inset 0 1px 0 0 rgba(255, 255, 255, 0.15); }
  .bp3-dark table.bp3-html-table.bp3-html-table-bordered tbody tr td:not(:first-child){
    -webkit-box-shadow:inset 1px 1px 0 0 rgba(255, 255, 255, 0.15);
            box-shadow:inset 1px 1px 0 0 rgba(255, 255, 255, 0.15); }

.bp3-dark table.bp3-html-table.bp3-html-table-bordered.bp3-html-table-striped tbody tr:not(:first-child) td{
  -webkit-box-shadow:inset 1px 0 0 0 rgba(255, 255, 255, 0.15);
          box-shadow:inset 1px 0 0 0 rgba(255, 255, 255, 0.15); }
  .bp3-dark table.bp3-html-table.bp3-html-table-bordered.bp3-html-table-striped tbody tr:not(:first-child) td:first-child{
    -webkit-box-shadow:none;
            box-shadow:none; }

.bp3-dark table.bp3-html-table.bp3-interactive tbody tr:hover td{
  background-color:rgba(92, 112, 128, 0.3);
  cursor:pointer; }

.bp3-dark table.bp3-html-table.bp3-interactive tbody tr:active td{
  background-color:rgba(92, 112, 128, 0.4); }

.bp3-key-combo{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center; }
  .bp3-key-combo > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-key-combo > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-key-combo::before,
  .bp3-key-combo > *{
    margin-right:5px; }
  .bp3-key-combo:empty::before,
  .bp3-key-combo > :last-child{
    margin-right:0; }

.bp3-hotkey-dialog{
  top:40px;
  padding-bottom:0; }
  .bp3-hotkey-dialog .bp3-dialog-body{
    margin:0;
    padding:0; }
  .bp3-hotkey-dialog .bp3-hotkey-label{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1; }

.bp3-hotkey-column{
  margin:auto;
  max-height:80vh;
  overflow-y:auto;
  padding:30px; }
  .bp3-hotkey-column .bp3-heading{
    margin-bottom:20px; }
    .bp3-hotkey-column .bp3-heading:not(:first-child){
      margin-top:40px; }

.bp3-hotkey{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:justify;
      -ms-flex-pack:justify;
          justify-content:space-between;
  margin-right:0;
  margin-left:0; }
  .bp3-hotkey:not(:last-child){
    margin-bottom:10px; }
.bp3-icon{
  display:inline-block;
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  vertical-align:text-bottom; }
  .bp3-icon:not(:empty)::before{
    content:"" !important;
    content:unset !important; }
  .bp3-icon > svg{
    display:block; }
    .bp3-icon > svg:not([fill]){
      fill:currentColor; }

.bp3-icon.bp3-intent-primary, .bp3-icon-standard.bp3-intent-primary, .bp3-icon-large.bp3-intent-primary{
  color:#106ba3; }
  .bp3-dark .bp3-icon.bp3-intent-primary, .bp3-dark .bp3-icon-standard.bp3-intent-primary, .bp3-dark .bp3-icon-large.bp3-intent-primary{
    color:#48aff0; }

.bp3-icon.bp3-intent-success, .bp3-icon-standard.bp3-intent-success, .bp3-icon-large.bp3-intent-success{
  color:#0d8050; }
  .bp3-dark .bp3-icon.bp3-intent-success, .bp3-dark .bp3-icon-standard.bp3-intent-success, .bp3-dark .bp3-icon-large.bp3-intent-success{
    color:#3dcc91; }

.bp3-icon.bp3-intent-warning, .bp3-icon-standard.bp3-intent-warning, .bp3-icon-large.bp3-intent-warning{
  color:#bf7326; }
  .bp3-dark .bp3-icon.bp3-intent-warning, .bp3-dark .bp3-icon-standard.bp3-intent-warning, .bp3-dark .bp3-icon-large.bp3-intent-warning{
    color:#ffb366; }

.bp3-icon.bp3-intent-danger, .bp3-icon-standard.bp3-intent-danger, .bp3-icon-large.bp3-intent-danger{
  color:#c23030; }
  .bp3-dark .bp3-icon.bp3-intent-danger, .bp3-dark .bp3-icon-standard.bp3-intent-danger, .bp3-dark .bp3-icon-large.bp3-intent-danger{
    color:#ff7373; }

span.bp3-icon-standard{
  line-height:1;
  font-family:"Icons16", sans-serif;
  font-size:16px;
  font-weight:400;
  font-style:normal;
  -moz-osx-font-smoothing:grayscale;
  -webkit-font-smoothing:antialiased;
  display:inline-block; }

span.bp3-icon-large{
  line-height:1;
  font-family:"Icons20", sans-serif;
  font-size:20px;
  font-weight:400;
  font-style:normal;
  -moz-osx-font-smoothing:grayscale;
  -webkit-font-smoothing:antialiased;
  display:inline-block; }

span.bp3-icon:empty{
  line-height:1;
  font-family:"Icons20";
  font-size:inherit;
  font-weight:400;
  font-style:normal; }
  span.bp3-icon:empty::before{
    -moz-osx-font-smoothing:grayscale;
    -webkit-font-smoothing:antialiased; }

.bp3-icon-add::before{
  content:""; }

.bp3-icon-add-column-left::before{
  content:""; }

.bp3-icon-add-column-right::before{
  content:""; }

.bp3-icon-add-row-bottom::before{
  content:""; }

.bp3-icon-add-row-top::before{
  content:""; }

.bp3-icon-add-to-artifact::before{
  content:""; }

.bp3-icon-add-to-folder::before{
  content:""; }

.bp3-icon-airplane::before{
  content:""; }

.bp3-icon-align-center::before{
  content:""; }

.bp3-icon-align-justify::before{
  content:""; }

.bp3-icon-align-left::before{
  content:""; }

.bp3-icon-align-right::before{
  content:""; }

.bp3-icon-alignment-bottom::before{
  content:""; }

.bp3-icon-alignment-horizontal-center::before{
  content:""; }

.bp3-icon-alignment-left::before{
  content:""; }

.bp3-icon-alignment-right::before{
  content:""; }

.bp3-icon-alignment-top::before{
  content:""; }

.bp3-icon-alignment-vertical-center::before{
  content:""; }

.bp3-icon-annotation::before{
  content:""; }

.bp3-icon-application::before{
  content:""; }

.bp3-icon-applications::before{
  content:""; }

.bp3-icon-archive::before{
  content:""; }

.bp3-icon-arrow-bottom-left::before{
  content:"↙"; }

.bp3-icon-arrow-bottom-right::before{
  content:"↘"; }

.bp3-icon-arrow-down::before{
  content:"↓"; }

.bp3-icon-arrow-left::before{
  content:"←"; }

.bp3-icon-arrow-right::before{
  content:"→"; }

.bp3-icon-arrow-top-left::before{
  content:"↖"; }

.bp3-icon-arrow-top-right::before{
  content:"↗"; }

.bp3-icon-arrow-up::before{
  content:"↑"; }

.bp3-icon-arrows-horizontal::before{
  content:"↔"; }

.bp3-icon-arrows-vertical::before{
  content:"↕"; }

.bp3-icon-asterisk::before{
  content:"*"; }

.bp3-icon-automatic-updates::before{
  content:""; }

.bp3-icon-badge::before{
  content:""; }

.bp3-icon-ban-circle::before{
  content:""; }

.bp3-icon-bank-account::before{
  content:""; }

.bp3-icon-barcode::before{
  content:""; }

.bp3-icon-blank::before{
  content:""; }

.bp3-icon-blocked-person::before{
  content:""; }

.bp3-icon-bold::before{
  content:""; }

.bp3-icon-book::before{
  content:""; }

.bp3-icon-bookmark::before{
  content:""; }

.bp3-icon-box::before{
  content:""; }

.bp3-icon-briefcase::before{
  content:""; }

.bp3-icon-bring-data::before{
  content:""; }

.bp3-icon-build::before{
  content:""; }

.bp3-icon-calculator::before{
  content:""; }

.bp3-icon-calendar::before{
  content:""; }

.bp3-icon-camera::before{
  content:""; }

.bp3-icon-caret-down::before{
  content:"⌄"; }

.bp3-icon-caret-left::before{
  content:"〈"; }

.bp3-icon-caret-right::before{
  content:"〉"; }

.bp3-icon-caret-up::before{
  content:"⌃"; }

.bp3-icon-cell-tower::before{
  content:""; }

.bp3-icon-changes::before{
  content:""; }

.bp3-icon-chart::before{
  content:""; }

.bp3-icon-chat::before{
  content:""; }

.bp3-icon-chevron-backward::before{
  content:""; }

.bp3-icon-chevron-down::before{
  content:""; }

.bp3-icon-chevron-forward::before{
  content:""; }

.bp3-icon-chevron-left::before{
  content:""; }

.bp3-icon-chevron-right::before{
  content:""; }

.bp3-icon-chevron-up::before{
  content:""; }

.bp3-icon-circle::before{
  content:""; }

.bp3-icon-circle-arrow-down::before{
  content:""; }

.bp3-icon-circle-arrow-left::before{
  content:""; }

.bp3-icon-circle-arrow-right::before{
  content:""; }

.bp3-icon-circle-arrow-up::before{
  content:""; }

.bp3-icon-citation::before{
  content:""; }

.bp3-icon-clean::before{
  content:""; }

.bp3-icon-clipboard::before{
  content:""; }

.bp3-icon-cloud::before{
  content:"☁"; }

.bp3-icon-cloud-download::before{
  content:""; }

.bp3-icon-cloud-upload::before{
  content:""; }

.bp3-icon-code::before{
  content:""; }

.bp3-icon-code-block::before{
  content:""; }

.bp3-icon-cog::before{
  content:""; }

.bp3-icon-collapse-all::before{
  content:""; }

.bp3-icon-column-layout::before{
  content:""; }

.bp3-icon-comment::before{
  content:""; }

.bp3-icon-comparison::before{
  content:""; }

.bp3-icon-compass::before{
  content:""; }

.bp3-icon-compressed::before{
  content:""; }

.bp3-icon-confirm::before{
  content:""; }

.bp3-icon-console::before{
  content:""; }

.bp3-icon-contrast::before{
  content:""; }

.bp3-icon-control::before{
  content:""; }

.bp3-icon-credit-card::before{
  content:""; }

.bp3-icon-cross::before{
  content:"✗"; }

.bp3-icon-crown::before{
  content:""; }

.bp3-icon-cube::before{
  content:""; }

.bp3-icon-cube-add::before{
  content:""; }

.bp3-icon-cube-remove::before{
  content:""; }

.bp3-icon-curved-range-chart::before{
  content:""; }

.bp3-icon-cut::before{
  content:""; }

.bp3-icon-dashboard::before{
  content:""; }

.bp3-icon-data-lineage::before{
  content:""; }

.bp3-icon-database::before{
  content:""; }

.bp3-icon-delete::before{
  content:""; }

.bp3-icon-delta::before{
  content:"Δ"; }

.bp3-icon-derive-column::before{
  content:""; }

.bp3-icon-desktop::before{
  content:""; }

.bp3-icon-diagram-tree::before{
  content:""; }

.bp3-icon-direction-left::before{
  content:""; }

.bp3-icon-direction-right::before{
  content:""; }

.bp3-icon-disable::before{
  content:""; }

.bp3-icon-document::before{
  content:""; }

.bp3-icon-document-open::before{
  content:""; }

.bp3-icon-document-share::before{
  content:""; }

.bp3-icon-dollar::before{
  content:"$"; }

.bp3-icon-dot::before{
  content:"•"; }

.bp3-icon-double-caret-horizontal::before{
  content:""; }

.bp3-icon-double-caret-vertical::before{
  content:""; }

.bp3-icon-double-chevron-down::before{
  content:""; }

.bp3-icon-double-chevron-left::before{
  content:""; }

.bp3-icon-double-chevron-right::before{
  content:""; }

.bp3-icon-double-chevron-up::before{
  content:""; }

.bp3-icon-doughnut-chart::before{
  content:""; }

.bp3-icon-download::before{
  content:""; }

.bp3-icon-drag-handle-horizontal::before{
  content:""; }

.bp3-icon-drag-handle-vertical::before{
  content:""; }

.bp3-icon-draw::before{
  content:""; }

.bp3-icon-drive-time::before{
  content:""; }

.bp3-icon-duplicate::before{
  content:""; }

.bp3-icon-edit::before{
  content:"✎"; }

.bp3-icon-eject::before{
  content:"⏏"; }

.bp3-icon-endorsed::before{
  content:""; }

.bp3-icon-envelope::before{
  content:"✉"; }

.bp3-icon-equals::before{
  content:""; }

.bp3-icon-eraser::before{
  content:""; }

.bp3-icon-error::before{
  content:""; }

.bp3-icon-euro::before{
  content:"€"; }

.bp3-icon-exchange::before{
  content:""; }

.bp3-icon-exclude-row::before{
  content:""; }

.bp3-icon-expand-all::before{
  content:""; }

.bp3-icon-export::before{
  content:""; }

.bp3-icon-eye-off::before{
  content:""; }

.bp3-icon-eye-on::before{
  content:""; }

.bp3-icon-eye-open::before{
  content:""; }

.bp3-icon-fast-backward::before{
  content:""; }

.bp3-icon-fast-forward::before{
  content:""; }

.bp3-icon-feed::before{
  content:""; }

.bp3-icon-feed-subscribed::before{
  content:""; }

.bp3-icon-film::before{
  content:""; }

.bp3-icon-filter::before{
  content:""; }

.bp3-icon-filter-keep::before{
  content:""; }

.bp3-icon-filter-list::before{
  content:""; }

.bp3-icon-filter-open::before{
  content:""; }

.bp3-icon-filter-remove::before{
  content:""; }

.bp3-icon-flag::before{
  content:"⚑"; }

.bp3-icon-flame::before{
  content:""; }

.bp3-icon-flash::before{
  content:""; }

.bp3-icon-floppy-disk::before{
  content:""; }

.bp3-icon-flow-branch::before{
  content:""; }

.bp3-icon-flow-end::before{
  content:""; }

.bp3-icon-flow-linear::before{
  content:""; }

.bp3-icon-flow-review::before{
  content:""; }

.bp3-icon-flow-review-branch::before{
  content:""; }

.bp3-icon-flows::before{
  content:""; }

.bp3-icon-folder-close::before{
  content:""; }

.bp3-icon-folder-new::before{
  content:""; }

.bp3-icon-folder-open::before{
  content:""; }

.bp3-icon-folder-shared::before{
  content:""; }

.bp3-icon-folder-shared-open::before{
  content:""; }

.bp3-icon-follower::before{
  content:""; }

.bp3-icon-following::before{
  content:""; }

.bp3-icon-font::before{
  content:""; }

.bp3-icon-fork::before{
  content:""; }

.bp3-icon-form::before{
  content:""; }

.bp3-icon-full-circle::before{
  content:""; }

.bp3-icon-full-stacked-chart::before{
  content:""; }

.bp3-icon-fullscreen::before{
  content:""; }

.bp3-icon-function::before{
  content:""; }

.bp3-icon-gantt-chart::before{
  content:""; }

.bp3-icon-geolocation::before{
  content:""; }

.bp3-icon-geosearch::before{
  content:""; }

.bp3-icon-git-branch::before{
  content:""; }

.bp3-icon-git-commit::before{
  content:""; }

.bp3-icon-git-merge::before{
  content:""; }

.bp3-icon-git-new-branch::before{
  content:""; }

.bp3-icon-git-pull::before{
  content:""; }

.bp3-icon-git-push::before{
  content:""; }

.bp3-icon-git-repo::before{
  content:""; }

.bp3-icon-glass::before{
  content:""; }

.bp3-icon-globe::before{
  content:""; }

.bp3-icon-globe-network::before{
  content:""; }

.bp3-icon-graph::before{
  content:""; }

.bp3-icon-graph-remove::before{
  content:""; }

.bp3-icon-greater-than::before{
  content:""; }

.bp3-icon-greater-than-or-equal-to::before{
  content:""; }

.bp3-icon-grid::before{
  content:""; }

.bp3-icon-grid-view::before{
  content:""; }

.bp3-icon-group-objects::before{
  content:""; }

.bp3-icon-grouped-bar-chart::before{
  content:""; }

.bp3-icon-hand::before{
  content:""; }

.bp3-icon-hand-down::before{
  content:""; }

.bp3-icon-hand-left::before{
  content:""; }

.bp3-icon-hand-right::before{
  content:""; }

.bp3-icon-hand-up::before{
  content:""; }

.bp3-icon-header::before{
  content:""; }

.bp3-icon-header-one::before{
  content:""; }

.bp3-icon-header-two::before{
  content:""; }

.bp3-icon-headset::before{
  content:""; }

.bp3-icon-heart::before{
  content:"♥"; }

.bp3-icon-heart-broken::before{
  content:""; }

.bp3-icon-heat-grid::before{
  content:""; }

.bp3-icon-heatmap::before{
  content:""; }

.bp3-icon-help::before{
  content:"?"; }

.bp3-icon-helper-management::before{
  content:""; }

.bp3-icon-highlight::before{
  content:""; }

.bp3-icon-history::before{
  content:""; }

.bp3-icon-home::before{
  content:"⌂"; }

.bp3-icon-horizontal-bar-chart::before{
  content:""; }

.bp3-icon-horizontal-bar-chart-asc::before{
  content:""; }

.bp3-icon-horizontal-bar-chart-desc::before{
  content:""; }

.bp3-icon-horizontal-distribution::before{
  content:""; }

.bp3-icon-id-number::before{
  content:""; }

.bp3-icon-image-rotate-left::before{
  content:""; }

.bp3-icon-image-rotate-right::before{
  content:""; }

.bp3-icon-import::before{
  content:""; }

.bp3-icon-inbox::before{
  content:""; }

.bp3-icon-inbox-filtered::before{
  content:""; }

.bp3-icon-inbox-geo::before{
  content:""; }

.bp3-icon-inbox-search::before{
  content:""; }

.bp3-icon-inbox-update::before{
  content:""; }

.bp3-icon-info-sign::before{
  content:"ℹ"; }

.bp3-icon-inheritance::before{
  content:""; }

.bp3-icon-inner-join::before{
  content:""; }

.bp3-icon-insert::before{
  content:""; }

.bp3-icon-intersection::before{
  content:""; }

.bp3-icon-ip-address::before{
  content:""; }

.bp3-icon-issue::before{
  content:""; }

.bp3-icon-issue-closed::before{
  content:""; }

.bp3-icon-issue-new::before{
  content:""; }

.bp3-icon-italic::before{
  content:""; }

.bp3-icon-join-table::before{
  content:""; }

.bp3-icon-key::before{
  content:""; }

.bp3-icon-key-backspace::before{
  content:""; }

.bp3-icon-key-command::before{
  content:""; }

.bp3-icon-key-control::before{
  content:""; }

.bp3-icon-key-delete::before{
  content:""; }

.bp3-icon-key-enter::before{
  content:""; }

.bp3-icon-key-escape::before{
  content:""; }

.bp3-icon-key-option::before{
  content:""; }

.bp3-icon-key-shift::before{
  content:""; }

.bp3-icon-key-tab::before{
  content:""; }

.bp3-icon-known-vehicle::before{
  content:""; }

.bp3-icon-label::before{
  content:""; }

.bp3-icon-layer::before{
  content:""; }

.bp3-icon-layers::before{
  content:""; }

.bp3-icon-layout::before{
  content:""; }

.bp3-icon-layout-auto::before{
  content:""; }

.bp3-icon-layout-balloon::before{
  content:""; }

.bp3-icon-layout-circle::before{
  content:""; }

.bp3-icon-layout-grid::before{
  content:""; }

.bp3-icon-layout-group-by::before{
  content:""; }

.bp3-icon-layout-hierarchy::before{
  content:""; }

.bp3-icon-layout-linear::before{
  content:""; }

.bp3-icon-layout-skew-grid::before{
  content:""; }

.bp3-icon-layout-sorted-clusters::before{
  content:""; }

.bp3-icon-learning::before{
  content:""; }

.bp3-icon-left-join::before{
  content:""; }

.bp3-icon-less-than::before{
  content:""; }

.bp3-icon-less-than-or-equal-to::before{
  content:""; }

.bp3-icon-lifesaver::before{
  content:""; }

.bp3-icon-lightbulb::before{
  content:""; }

.bp3-icon-link::before{
  content:""; }

.bp3-icon-list::before{
  content:"☰"; }

.bp3-icon-list-columns::before{
  content:""; }

.bp3-icon-list-detail-view::before{
  content:""; }

.bp3-icon-locate::before{
  content:""; }

.bp3-icon-lock::before{
  content:""; }

.bp3-icon-log-in::before{
  content:""; }

.bp3-icon-log-out::before{
  content:""; }

.bp3-icon-manual::before{
  content:""; }

.bp3-icon-manually-entered-data::before{
  content:""; }

.bp3-icon-map::before{
  content:""; }

.bp3-icon-map-create::before{
  content:""; }

.bp3-icon-map-marker::before{
  content:""; }

.bp3-icon-maximize::before{
  content:""; }

.bp3-icon-media::before{
  content:""; }

.bp3-icon-menu::before{
  content:""; }

.bp3-icon-menu-closed::before{
  content:""; }

.bp3-icon-menu-open::before{
  content:""; }

.bp3-icon-merge-columns::before{
  content:""; }

.bp3-icon-merge-links::before{
  content:""; }

.bp3-icon-minimize::before{
  content:""; }

.bp3-icon-minus::before{
  content:"−"; }

.bp3-icon-mobile-phone::before{
  content:""; }

.bp3-icon-mobile-video::before{
  content:""; }

.bp3-icon-moon::before{
  content:""; }

.bp3-icon-more::before{
  content:""; }

.bp3-icon-mountain::before{
  content:""; }

.bp3-icon-move::before{
  content:""; }

.bp3-icon-mugshot::before{
  content:""; }

.bp3-icon-multi-select::before{
  content:""; }

.bp3-icon-music::before{
  content:""; }

.bp3-icon-new-drawing::before{
  content:""; }

.bp3-icon-new-grid-item::before{
  content:""; }

.bp3-icon-new-layer::before{
  content:""; }

.bp3-icon-new-layers::before{
  content:""; }

.bp3-icon-new-link::before{
  content:""; }

.bp3-icon-new-object::before{
  content:""; }

.bp3-icon-new-person::before{
  content:""; }

.bp3-icon-new-prescription::before{
  content:""; }

.bp3-icon-new-text-box::before{
  content:""; }

.bp3-icon-ninja::before{
  content:""; }

.bp3-icon-not-equal-to::before{
  content:""; }

.bp3-icon-notifications::before{
  content:""; }

.bp3-icon-notifications-updated::before{
  content:""; }

.bp3-icon-numbered-list::before{
  content:""; }

.bp3-icon-numerical::before{
  content:""; }

.bp3-icon-office::before{
  content:""; }

.bp3-icon-offline::before{
  content:""; }

.bp3-icon-oil-field::before{
  content:""; }

.bp3-icon-one-column::before{
  content:""; }

.bp3-icon-outdated::before{
  content:""; }

.bp3-icon-page-layout::before{
  content:""; }

.bp3-icon-panel-stats::before{
  content:""; }

.bp3-icon-panel-table::before{
  content:""; }

.bp3-icon-paperclip::before{
  content:""; }

.bp3-icon-paragraph::before{
  content:""; }

.bp3-icon-path::before{
  content:""; }

.bp3-icon-path-search::before{
  content:""; }

.bp3-icon-pause::before{
  content:""; }

.bp3-icon-people::before{
  content:""; }

.bp3-icon-percentage::before{
  content:""; }

.bp3-icon-person::before{
  content:""; }

.bp3-icon-phone::before{
  content:"☎"; }

.bp3-icon-pie-chart::before{
  content:""; }

.bp3-icon-pin::before{
  content:""; }

.bp3-icon-pivot::before{
  content:""; }

.bp3-icon-pivot-table::before{
  content:""; }

.bp3-icon-play::before{
  content:""; }

.bp3-icon-plus::before{
  content:"+"; }

.bp3-icon-polygon-filter::before{
  content:""; }

.bp3-icon-power::before{
  content:""; }

.bp3-icon-predictive-analysis::before{
  content:""; }

.bp3-icon-prescription::before{
  content:""; }

.bp3-icon-presentation::before{
  content:""; }

.bp3-icon-print::before{
  content:"⎙"; }

.bp3-icon-projects::before{
  content:""; }

.bp3-icon-properties::before{
  content:""; }

.bp3-icon-property::before{
  content:""; }

.bp3-icon-publish-function::before{
  content:""; }

.bp3-icon-pulse::before{
  content:""; }

.bp3-icon-random::before{
  content:""; }

.bp3-icon-record::before{
  content:""; }

.bp3-icon-redo::before{
  content:""; }

.bp3-icon-refresh::before{
  content:""; }

.bp3-icon-regression-chart::before{
  content:""; }

.bp3-icon-remove::before{
  content:""; }

.bp3-icon-remove-column::before{
  content:""; }

.bp3-icon-remove-column-left::before{
  content:""; }

.bp3-icon-remove-column-right::before{
  content:""; }

.bp3-icon-remove-row-bottom::before{
  content:""; }

.bp3-icon-remove-row-top::before{
  content:""; }

.bp3-icon-repeat::before{
  content:""; }

.bp3-icon-reset::before{
  content:""; }

.bp3-icon-resolve::before{
  content:""; }

.bp3-icon-rig::before{
  content:""; }

.bp3-icon-right-join::before{
  content:""; }

.bp3-icon-ring::before{
  content:""; }

.bp3-icon-rotate-document::before{
  content:""; }

.bp3-icon-rotate-page::before{
  content:""; }

.bp3-icon-satellite::before{
  content:""; }

.bp3-icon-saved::before{
  content:""; }

.bp3-icon-scatter-plot::before{
  content:""; }

.bp3-icon-search::before{
  content:""; }

.bp3-icon-search-around::before{
  content:""; }

.bp3-icon-search-template::before{
  content:""; }

.bp3-icon-search-text::before{
  content:""; }

.bp3-icon-segmented-control::before{
  content:""; }

.bp3-icon-select::before{
  content:""; }

.bp3-icon-selection::before{
  content:"⦿"; }

.bp3-icon-send-to::before{
  content:""; }

.bp3-icon-send-to-graph::before{
  content:""; }

.bp3-icon-send-to-map::before{
  content:""; }

.bp3-icon-series-add::before{
  content:""; }

.bp3-icon-series-configuration::before{
  content:""; }

.bp3-icon-series-derived::before{
  content:""; }

.bp3-icon-series-filtered::before{
  content:""; }

.bp3-icon-series-search::before{
  content:""; }

.bp3-icon-settings::before{
  content:""; }

.bp3-icon-share::before{
  content:""; }

.bp3-icon-shield::before{
  content:""; }

.bp3-icon-shop::before{
  content:""; }

.bp3-icon-shopping-cart::before{
  content:""; }

.bp3-icon-signal-search::before{
  content:""; }

.bp3-icon-sim-card::before{
  content:""; }

.bp3-icon-slash::before{
  content:""; }

.bp3-icon-small-cross::before{
  content:""; }

.bp3-icon-small-minus::before{
  content:""; }

.bp3-icon-small-plus::before{
  content:""; }

.bp3-icon-small-tick::before{
  content:""; }

.bp3-icon-snowflake::before{
  content:""; }

.bp3-icon-social-media::before{
  content:""; }

.bp3-icon-sort::before{
  content:""; }

.bp3-icon-sort-alphabetical::before{
  content:""; }

.bp3-icon-sort-alphabetical-desc::before{
  content:""; }

.bp3-icon-sort-asc::before{
  content:""; }

.bp3-icon-sort-desc::before{
  content:""; }

.bp3-icon-sort-numerical::before{
  content:""; }

.bp3-icon-sort-numerical-desc::before{
  content:""; }

.bp3-icon-split-columns::before{
  content:""; }

.bp3-icon-square::before{
  content:""; }

.bp3-icon-stacked-chart::before{
  content:""; }

.bp3-icon-star::before{
  content:"★"; }

.bp3-icon-star-empty::before{
  content:"☆"; }

.bp3-icon-step-backward::before{
  content:""; }

.bp3-icon-step-chart::before{
  content:""; }

.bp3-icon-step-forward::before{
  content:""; }

.bp3-icon-stop::before{
  content:""; }

.bp3-icon-stopwatch::before{
  content:""; }

.bp3-icon-strikethrough::before{
  content:""; }

.bp3-icon-style::before{
  content:""; }

.bp3-icon-swap-horizontal::before{
  content:""; }

.bp3-icon-swap-vertical::before{
  content:""; }

.bp3-icon-symbol-circle::before{
  content:""; }

.bp3-icon-symbol-cross::before{
  content:""; }

.bp3-icon-symbol-diamond::before{
  content:""; }

.bp3-icon-symbol-square::before{
  content:""; }

.bp3-icon-symbol-triangle-down::before{
  content:""; }

.bp3-icon-symbol-triangle-up::before{
  content:""; }

.bp3-icon-tag::before{
  content:""; }

.bp3-icon-take-action::before{
  content:""; }

.bp3-icon-taxi::before{
  content:""; }

.bp3-icon-text-highlight::before{
  content:""; }

.bp3-icon-th::before{
  content:""; }

.bp3-icon-th-derived::before{
  content:""; }

.bp3-icon-th-disconnect::before{
  content:""; }

.bp3-icon-th-filtered::before{
  content:""; }

.bp3-icon-th-list::before{
  content:""; }

.bp3-icon-thumbs-down::before{
  content:""; }

.bp3-icon-thumbs-up::before{
  content:""; }

.bp3-icon-tick::before{
  content:"✓"; }

.bp3-icon-tick-circle::before{
  content:""; }

.bp3-icon-time::before{
  content:"⏲"; }

.bp3-icon-timeline-area-chart::before{
  content:""; }

.bp3-icon-timeline-bar-chart::before{
  content:""; }

.bp3-icon-timeline-events::before{
  content:""; }

.bp3-icon-timeline-line-chart::before{
  content:""; }

.bp3-icon-tint::before{
  content:""; }

.bp3-icon-torch::before{
  content:""; }

.bp3-icon-tractor::before{
  content:""; }

.bp3-icon-train::before{
  content:""; }

.bp3-icon-translate::before{
  content:""; }

.bp3-icon-trash::before{
  content:""; }

.bp3-icon-tree::before{
  content:""; }

.bp3-icon-trending-down::before{
  content:""; }

.bp3-icon-trending-up::before{
  content:""; }

.bp3-icon-truck::before{
  content:""; }

.bp3-icon-two-columns::before{
  content:""; }

.bp3-icon-unarchive::before{
  content:""; }

.bp3-icon-underline::before{
  content:"⎁"; }

.bp3-icon-undo::before{
  content:"⎌"; }

.bp3-icon-ungroup-objects::before{
  content:""; }

.bp3-icon-unknown-vehicle::before{
  content:""; }

.bp3-icon-unlock::before{
  content:""; }

.bp3-icon-unpin::before{
  content:""; }

.bp3-icon-unresolve::before{
  content:""; }

.bp3-icon-updated::before{
  content:""; }

.bp3-icon-upload::before{
  content:""; }

.bp3-icon-user::before{
  content:""; }

.bp3-icon-variable::before{
  content:""; }

.bp3-icon-vertical-bar-chart-asc::before{
  content:""; }

.bp3-icon-vertical-bar-chart-desc::before{
  content:""; }

.bp3-icon-vertical-distribution::before{
  content:""; }

.bp3-icon-video::before{
  content:""; }

.bp3-icon-volume-down::before{
  content:""; }

.bp3-icon-volume-off::before{
  content:""; }

.bp3-icon-volume-up::before{
  content:""; }

.bp3-icon-walk::before{
  content:""; }

.bp3-icon-warning-sign::before{
  content:""; }

.bp3-icon-waterfall-chart::before{
  content:""; }

.bp3-icon-widget::before{
  content:""; }

.bp3-icon-widget-button::before{
  content:""; }

.bp3-icon-widget-footer::before{
  content:""; }

.bp3-icon-widget-header::before{
  content:""; }

.bp3-icon-wrench::before{
  content:""; }

.bp3-icon-zoom-in::before{
  content:""; }

.bp3-icon-zoom-out::before{
  content:""; }

.bp3-icon-zoom-to-fit::before{
  content:""; }
.bp3-submenu > .bp3-popover-wrapper{
  display:block; }

.bp3-submenu .bp3-popover-target{
  display:block; }

.bp3-submenu.bp3-popover{
  -webkit-box-shadow:none;
          box-shadow:none;
  padding:0 5px; }
  .bp3-submenu.bp3-popover > .bp3-popover-content{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2); }
  .bp3-dark .bp3-submenu.bp3-popover, .bp3-submenu.bp3-popover.bp3-dark{
    -webkit-box-shadow:none;
            box-shadow:none; }
    .bp3-dark .bp3-submenu.bp3-popover > .bp3-popover-content, .bp3-submenu.bp3-popover.bp3-dark > .bp3-popover-content{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4); }
.bp3-menu{
  margin:0;
  border-radius:3px;
  background:#ffffff;
  min-width:180px;
  padding:5px;
  list-style:none;
  text-align:left;
  color:#182026; }

.bp3-menu-divider{
  display:block;
  margin:5px;
  border-top:1px solid rgba(16, 22, 26, 0.15); }
  .bp3-dark .bp3-menu-divider{
    border-top-color:rgba(255, 255, 255, 0.15); }

.bp3-menu-item{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:start;
      -ms-flex-align:start;
          align-items:flex-start;
  border-radius:2px;
  padding:5px 7px;
  text-decoration:none;
  line-height:20px;
  color:inherit;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-menu-item > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-menu-item > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-menu-item::before,
  .bp3-menu-item > *{
    margin-right:7px; }
  .bp3-menu-item:empty::before,
  .bp3-menu-item > :last-child{
    margin-right:0; }
  .bp3-menu-item > .bp3-fill{
    word-break:break-word; }
  .bp3-menu-item:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-menu-item{
    background-color:rgba(167, 182, 194, 0.3);
    cursor:pointer;
    text-decoration:none; }
  .bp3-menu-item.bp3-disabled{
    background-color:inherit;
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-dark .bp3-menu-item{
    color:inherit; }
    .bp3-dark .bp3-menu-item:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-menu-item{
      background-color:rgba(138, 155, 168, 0.15);
      color:inherit; }
    .bp3-dark .bp3-menu-item.bp3-disabled{
      background-color:inherit;
      color:rgba(167, 182, 194, 0.6); }
  .bp3-menu-item.bp3-intent-primary{
    color:#106ba3; }
    .bp3-menu-item.bp3-intent-primary .bp3-icon{
      color:inherit; }
    .bp3-menu-item.bp3-intent-primary::before, .bp3-menu-item.bp3-intent-primary::after,
    .bp3-menu-item.bp3-intent-primary .bp3-menu-item-label{
      color:#106ba3; }
    .bp3-menu-item.bp3-intent-primary:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item, .bp3-menu-item.bp3-intent-primary.bp3-active{
      background-color:#137cbd; }
    .bp3-menu-item.bp3-intent-primary:active{
      background-color:#106ba3; }
    .bp3-menu-item.bp3-intent-primary:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item, .bp3-menu-item.bp3-intent-primary:hover::before, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item::before, .bp3-menu-item.bp3-intent-primary:hover::after, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item::after,
    .bp3-menu-item.bp3-intent-primary:hover .bp3-menu-item-label,
    .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item .bp3-menu-item-label, .bp3-menu-item.bp3-intent-primary:active, .bp3-menu-item.bp3-intent-primary:active::before, .bp3-menu-item.bp3-intent-primary:active::after,
    .bp3-menu-item.bp3-intent-primary:active .bp3-menu-item-label, .bp3-menu-item.bp3-intent-primary.bp3-active, .bp3-menu-item.bp3-intent-primary.bp3-active::before, .bp3-menu-item.bp3-intent-primary.bp3-active::after,
    .bp3-menu-item.bp3-intent-primary.bp3-active .bp3-menu-item-label{
      color:#ffffff; }
  .bp3-menu-item.bp3-intent-success{
    color:#0d8050; }
    .bp3-menu-item.bp3-intent-success .bp3-icon{
      color:inherit; }
    .bp3-menu-item.bp3-intent-success::before, .bp3-menu-item.bp3-intent-success::after,
    .bp3-menu-item.bp3-intent-success .bp3-menu-item-label{
      color:#0d8050; }
    .bp3-menu-item.bp3-intent-success:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item, .bp3-menu-item.bp3-intent-success.bp3-active{
      background-color:#0f9960; }
    .bp3-menu-item.bp3-intent-success:active{
      background-color:#0d8050; }
    .bp3-menu-item.bp3-intent-success:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item, .bp3-menu-item.bp3-intent-success:hover::before, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item::before, .bp3-menu-item.bp3-intent-success:hover::after, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item::after,
    .bp3-menu-item.bp3-intent-success:hover .bp3-menu-item-label,
    .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item .bp3-menu-item-label, .bp3-menu-item.bp3-intent-success:active, .bp3-menu-item.bp3-intent-success:active::before, .bp3-menu-item.bp3-intent-success:active::after,
    .bp3-menu-item.bp3-intent-success:active .bp3-menu-item-label, .bp3-menu-item.bp3-intent-success.bp3-active, .bp3-menu-item.bp3-intent-success.bp3-active::before, .bp3-menu-item.bp3-intent-success.bp3-active::after,
    .bp3-menu-item.bp3-intent-success.bp3-active .bp3-menu-item-label{
      color:#ffffff; }
  .bp3-menu-item.bp3-intent-warning{
    color:#bf7326; }
    .bp3-menu-item.bp3-intent-warning .bp3-icon{
      color:inherit; }
    .bp3-menu-item.bp3-intent-warning::before, .bp3-menu-item.bp3-intent-warning::after,
    .bp3-menu-item.bp3-intent-warning .bp3-menu-item-label{
      color:#bf7326; }
    .bp3-menu-item.bp3-intent-warning:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item, .bp3-menu-item.bp3-intent-warning.bp3-active{
      background-color:#d9822b; }
    .bp3-menu-item.bp3-intent-warning:active{
      background-color:#bf7326; }
    .bp3-menu-item.bp3-intent-warning:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item, .bp3-menu-item.bp3-intent-warning:hover::before, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item::before, .bp3-menu-item.bp3-intent-warning:hover::after, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item::after,
    .bp3-menu-item.bp3-intent-warning:hover .bp3-menu-item-label,
    .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item .bp3-menu-item-label, .bp3-menu-item.bp3-intent-warning:active, .bp3-menu-item.bp3-intent-warning:active::before, .bp3-menu-item.bp3-intent-warning:active::after,
    .bp3-menu-item.bp3-intent-warning:active .bp3-menu-item-label, .bp3-menu-item.bp3-intent-warning.bp3-active, .bp3-menu-item.bp3-intent-warning.bp3-active::before, .bp3-menu-item.bp3-intent-warning.bp3-active::after,
    .bp3-menu-item.bp3-intent-warning.bp3-active .bp3-menu-item-label{
      color:#ffffff; }
  .bp3-menu-item.bp3-intent-danger{
    color:#c23030; }
    .bp3-menu-item.bp3-intent-danger .bp3-icon{
      color:inherit; }
    .bp3-menu-item.bp3-intent-danger::before, .bp3-menu-item.bp3-intent-danger::after,
    .bp3-menu-item.bp3-intent-danger .bp3-menu-item-label{
      color:#c23030; }
    .bp3-menu-item.bp3-intent-danger:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item, .bp3-menu-item.bp3-intent-danger.bp3-active{
      background-color:#db3737; }
    .bp3-menu-item.bp3-intent-danger:active{
      background-color:#c23030; }
    .bp3-menu-item.bp3-intent-danger:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item, .bp3-menu-item.bp3-intent-danger:hover::before, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item::before, .bp3-menu-item.bp3-intent-danger:hover::after, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item::after,
    .bp3-menu-item.bp3-intent-danger:hover .bp3-menu-item-label,
    .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item .bp3-menu-item-label, .bp3-menu-item.bp3-intent-danger:active, .bp3-menu-item.bp3-intent-danger:active::before, .bp3-menu-item.bp3-intent-danger:active::after,
    .bp3-menu-item.bp3-intent-danger:active .bp3-menu-item-label, .bp3-menu-item.bp3-intent-danger.bp3-active, .bp3-menu-item.bp3-intent-danger.bp3-active::before, .bp3-menu-item.bp3-intent-danger.bp3-active::after,
    .bp3-menu-item.bp3-intent-danger.bp3-active .bp3-menu-item-label{
      color:#ffffff; }
  .bp3-menu-item::before{
    line-height:1;
    font-family:"Icons16", sans-serif;
    font-size:16px;
    font-weight:400;
    font-style:normal;
    -moz-osx-font-smoothing:grayscale;
    -webkit-font-smoothing:antialiased;
    margin-right:7px; }
  .bp3-menu-item::before,
  .bp3-menu-item > .bp3-icon{
    margin-top:2px;
    color:#5c7080; }
  .bp3-menu-item .bp3-menu-item-label{
    color:#5c7080; }
  .bp3-menu-item:hover, .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-menu-item{
    color:inherit; }
  .bp3-menu-item.bp3-active, .bp3-menu-item:active{
    background-color:rgba(115, 134, 148, 0.3); }
  .bp3-menu-item.bp3-disabled{
    outline:none !important;
    background-color:inherit !important;
    cursor:not-allowed !important;
    color:rgba(92, 112, 128, 0.6) !important; }
    .bp3-menu-item.bp3-disabled::before,
    .bp3-menu-item.bp3-disabled > .bp3-icon,
    .bp3-menu-item.bp3-disabled .bp3-menu-item-label{
      color:rgba(92, 112, 128, 0.6) !important; }
  .bp3-large .bp3-menu-item{
    padding:9px 7px;
    line-height:22px;
    font-size:16px; }
    .bp3-large .bp3-menu-item .bp3-icon{
      margin-top:3px; }
    .bp3-large .bp3-menu-item::before{
      line-height:1;
      font-family:"Icons20", sans-serif;
      font-size:20px;
      font-weight:400;
      font-style:normal;
      -moz-osx-font-smoothing:grayscale;
      -webkit-font-smoothing:antialiased;
      margin-top:1px;
      margin-right:10px; }

button.bp3-menu-item{
  border:none;
  background:none;
  width:100%;
  text-align:left; }
.bp3-menu-header{
  display:block;
  margin:5px;
  border-top:1px solid rgba(16, 22, 26, 0.15);
  cursor:default;
  padding-left:2px; }
  .bp3-dark .bp3-menu-header{
    border-top-color:rgba(255, 255, 255, 0.15); }
  .bp3-menu-header:first-of-type{
    border-top:none; }
  .bp3-menu-header > h6{
    color:#182026;
    font-weight:600;
    overflow:hidden;
    text-overflow:ellipsis;
    white-space:nowrap;
    word-wrap:normal;
    margin:0;
    padding:10px 7px 0 1px;
    line-height:17px; }
    .bp3-dark .bp3-menu-header > h6{
      color:#f5f8fa; }
  .bp3-menu-header:first-of-type > h6{
    padding-top:0; }
  .bp3-large .bp3-menu-header > h6{
    padding-top:15px;
    padding-bottom:5px;
    font-size:18px; }
  .bp3-large .bp3-menu-header:first-of-type > h6{
    padding-top:0; }

.bp3-dark .bp3-menu{
  background:#30404d;
  color:#f5f8fa; }

.bp3-dark .bp3-menu-item.bp3-intent-primary{
  color:#48aff0; }
  .bp3-dark .bp3-menu-item.bp3-intent-primary .bp3-icon{
    color:inherit; }
  .bp3-dark .bp3-menu-item.bp3-intent-primary::before, .bp3-dark .bp3-menu-item.bp3-intent-primary::after,
  .bp3-dark .bp3-menu-item.bp3-intent-primary .bp3-menu-item-label{
    color:#48aff0; }
  .bp3-dark .bp3-menu-item.bp3-intent-primary:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-primary.bp3-active{
    background-color:#137cbd; }
  .bp3-dark .bp3-menu-item.bp3-intent-primary:active{
    background-color:#106ba3; }
  .bp3-dark .bp3-menu-item.bp3-intent-primary:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-primary:hover::before, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item::before, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item::before, .bp3-dark .bp3-menu-item.bp3-intent-primary:hover::after, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item::after, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item::after,
  .bp3-dark .bp3-menu-item.bp3-intent-primary:hover .bp3-menu-item-label,
  .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item .bp3-menu-item-label,
  .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-primary.bp3-menu-item .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-primary:active, .bp3-dark .bp3-menu-item.bp3-intent-primary:active::before, .bp3-dark .bp3-menu-item.bp3-intent-primary:active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-primary:active .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-primary.bp3-active, .bp3-dark .bp3-menu-item.bp3-intent-primary.bp3-active::before, .bp3-dark .bp3-menu-item.bp3-intent-primary.bp3-active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-primary.bp3-active .bp3-menu-item-label{
    color:#ffffff; }

.bp3-dark .bp3-menu-item.bp3-intent-success{
  color:#3dcc91; }
  .bp3-dark .bp3-menu-item.bp3-intent-success .bp3-icon{
    color:inherit; }
  .bp3-dark .bp3-menu-item.bp3-intent-success::before, .bp3-dark .bp3-menu-item.bp3-intent-success::after,
  .bp3-dark .bp3-menu-item.bp3-intent-success .bp3-menu-item-label{
    color:#3dcc91; }
  .bp3-dark .bp3-menu-item.bp3-intent-success:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-success.bp3-active{
    background-color:#0f9960; }
  .bp3-dark .bp3-menu-item.bp3-intent-success:active{
    background-color:#0d8050; }
  .bp3-dark .bp3-menu-item.bp3-intent-success:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-success:hover::before, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item::before, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item::before, .bp3-dark .bp3-menu-item.bp3-intent-success:hover::after, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item::after, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item::after,
  .bp3-dark .bp3-menu-item.bp3-intent-success:hover .bp3-menu-item-label,
  .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item .bp3-menu-item-label,
  .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-success.bp3-menu-item .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-success:active, .bp3-dark .bp3-menu-item.bp3-intent-success:active::before, .bp3-dark .bp3-menu-item.bp3-intent-success:active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-success:active .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-success.bp3-active, .bp3-dark .bp3-menu-item.bp3-intent-success.bp3-active::before, .bp3-dark .bp3-menu-item.bp3-intent-success.bp3-active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-success.bp3-active .bp3-menu-item-label{
    color:#ffffff; }

.bp3-dark .bp3-menu-item.bp3-intent-warning{
  color:#ffb366; }
  .bp3-dark .bp3-menu-item.bp3-intent-warning .bp3-icon{
    color:inherit; }
  .bp3-dark .bp3-menu-item.bp3-intent-warning::before, .bp3-dark .bp3-menu-item.bp3-intent-warning::after,
  .bp3-dark .bp3-menu-item.bp3-intent-warning .bp3-menu-item-label{
    color:#ffb366; }
  .bp3-dark .bp3-menu-item.bp3-intent-warning:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-warning.bp3-active{
    background-color:#d9822b; }
  .bp3-dark .bp3-menu-item.bp3-intent-warning:active{
    background-color:#bf7326; }
  .bp3-dark .bp3-menu-item.bp3-intent-warning:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-warning:hover::before, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item::before, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item::before, .bp3-dark .bp3-menu-item.bp3-intent-warning:hover::after, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item::after, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item::after,
  .bp3-dark .bp3-menu-item.bp3-intent-warning:hover .bp3-menu-item-label,
  .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item .bp3-menu-item-label,
  .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-warning.bp3-menu-item .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-warning:active, .bp3-dark .bp3-menu-item.bp3-intent-warning:active::before, .bp3-dark .bp3-menu-item.bp3-intent-warning:active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-warning:active .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-warning.bp3-active, .bp3-dark .bp3-menu-item.bp3-intent-warning.bp3-active::before, .bp3-dark .bp3-menu-item.bp3-intent-warning.bp3-active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-warning.bp3-active .bp3-menu-item-label{
    color:#ffffff; }

.bp3-dark .bp3-menu-item.bp3-intent-danger{
  color:#ff7373; }
  .bp3-dark .bp3-menu-item.bp3-intent-danger .bp3-icon{
    color:inherit; }
  .bp3-dark .bp3-menu-item.bp3-intent-danger::before, .bp3-dark .bp3-menu-item.bp3-intent-danger::after,
  .bp3-dark .bp3-menu-item.bp3-intent-danger .bp3-menu-item-label{
    color:#ff7373; }
  .bp3-dark .bp3-menu-item.bp3-intent-danger:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-danger.bp3-active{
    background-color:#db3737; }
  .bp3-dark .bp3-menu-item.bp3-intent-danger:active{
    background-color:#c23030; }
  .bp3-dark .bp3-menu-item.bp3-intent-danger:hover, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item, .bp3-dark .bp3-menu-item.bp3-intent-danger:hover::before, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item::before, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item::before, .bp3-dark .bp3-menu-item.bp3-intent-danger:hover::after, .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item::after, .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item::after,
  .bp3-dark .bp3-menu-item.bp3-intent-danger:hover .bp3-menu-item-label,
  .bp3-dark .bp3-submenu .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item .bp3-menu-item-label,
  .bp3-submenu .bp3-dark .bp3-popover-target.bp3-popover-open > .bp3-intent-danger.bp3-menu-item .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-danger:active, .bp3-dark .bp3-menu-item.bp3-intent-danger:active::before, .bp3-dark .bp3-menu-item.bp3-intent-danger:active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-danger:active .bp3-menu-item-label, .bp3-dark .bp3-menu-item.bp3-intent-danger.bp3-active, .bp3-dark .bp3-menu-item.bp3-intent-danger.bp3-active::before, .bp3-dark .bp3-menu-item.bp3-intent-danger.bp3-active::after,
  .bp3-dark .bp3-menu-item.bp3-intent-danger.bp3-active .bp3-menu-item-label{
    color:#ffffff; }

.bp3-dark .bp3-menu-item::before,
.bp3-dark .bp3-menu-item > .bp3-icon{
  color:#a7b6c2; }

.bp3-dark .bp3-menu-item .bp3-menu-item-label{
  color:#a7b6c2; }

.bp3-dark .bp3-menu-item.bp3-active, .bp3-dark .bp3-menu-item:active{
  background-color:rgba(138, 155, 168, 0.3); }

.bp3-dark .bp3-menu-item.bp3-disabled{
  color:rgba(167, 182, 194, 0.6) !important; }
  .bp3-dark .bp3-menu-item.bp3-disabled::before,
  .bp3-dark .bp3-menu-item.bp3-disabled > .bp3-icon,
  .bp3-dark .bp3-menu-item.bp3-disabled .bp3-menu-item-label{
    color:rgba(167, 182, 194, 0.6) !important; }

.bp3-dark .bp3-menu-divider,
.bp3-dark .bp3-menu-header{
  border-color:rgba(255, 255, 255, 0.15); }

.bp3-dark .bp3-menu-header > h6{
  color:#f5f8fa; }

.bp3-label .bp3-menu{
  margin-top:5px; }
.bp3-navbar{
  position:relative;
  z-index:10;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.2);
  background-color:#ffffff;
  width:100%;
  height:50px;
  padding:0 15px; }
  .bp3-navbar.bp3-dark,
  .bp3-dark .bp3-navbar{
    background-color:#394b59; }
  .bp3-navbar.bp3-dark{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4); }
  .bp3-dark .bp3-navbar{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 0 0 rgba(16, 22, 26, 0), 0 1px 1px rgba(16, 22, 26, 0.4); }
  .bp3-navbar.bp3-fixed-top{
    position:fixed;
    top:0;
    right:0;
    left:0; }

.bp3-navbar-heading{
  margin-right:15px;
  font-size:16px; }

.bp3-navbar-group{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  height:50px; }
  .bp3-navbar-group.bp3-align-left{
    float:left; }
  .bp3-navbar-group.bp3-align-right{
    float:right; }

.bp3-navbar-divider{
  margin:0 10px;
  border-left:1px solid rgba(16, 22, 26, 0.15);
  height:20px; }
  .bp3-dark .bp3-navbar-divider{
    border-left-color:rgba(255, 255, 255, 0.15); }
.bp3-non-ideal-state{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:vertical;
  -webkit-box-direction:normal;
      -ms-flex-direction:column;
          flex-direction:column;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:center;
      -ms-flex-pack:center;
          justify-content:center;
  width:100%;
  height:100%;
  text-align:center; }
  .bp3-non-ideal-state > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-non-ideal-state > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-non-ideal-state::before,
  .bp3-non-ideal-state > *{
    margin-bottom:20px; }
  .bp3-non-ideal-state:empty::before,
  .bp3-non-ideal-state > :last-child{
    margin-bottom:0; }
  .bp3-non-ideal-state > *{
    max-width:400px; }

.bp3-non-ideal-state-visual{
  color:rgba(92, 112, 128, 0.6);
  font-size:60px; }
  .bp3-dark .bp3-non-ideal-state-visual{
    color:rgba(167, 182, 194, 0.6); }

.bp3-overflow-list{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -ms-flex-wrap:nowrap;
      flex-wrap:nowrap;
  min-width:0; }

.bp3-overflow-list-spacer{
  -ms-flex-negative:1;
      flex-shrink:1;
  width:1px; }

body.bp3-overlay-open{
  overflow:hidden; }

.bp3-overlay{
  position:static;
  top:0;
  right:0;
  bottom:0;
  left:0;
  z-index:20; }
  .bp3-overlay:not(.bp3-overlay-open){
    pointer-events:none; }
  .bp3-overlay.bp3-overlay-container{
    position:fixed;
    overflow:hidden; }
    .bp3-overlay.bp3-overlay-container.bp3-overlay-inline{
      position:absolute; }
  .bp3-overlay.bp3-overlay-scroll-container{
    position:fixed;
    overflow:auto; }
    .bp3-overlay.bp3-overlay-scroll-container.bp3-overlay-inline{
      position:absolute; }
  .bp3-overlay.bp3-overlay-inline{
    display:inline;
    overflow:visible; }

.bp3-overlay-content{
  position:fixed;
  z-index:20; }
  .bp3-overlay-inline .bp3-overlay-content,
  .bp3-overlay-scroll-container .bp3-overlay-content{
    position:absolute; }

.bp3-overlay-backdrop{
  position:fixed;
  top:0;
  right:0;
  bottom:0;
  left:0;
  opacity:1;
  z-index:20;
  background-color:rgba(16, 22, 26, 0.7);
  overflow:auto;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-overlay-backdrop.bp3-overlay-enter, .bp3-overlay-backdrop.bp3-overlay-appear{
    opacity:0; }
  .bp3-overlay-backdrop.bp3-overlay-enter-active, .bp3-overlay-backdrop.bp3-overlay-appear-active{
    opacity:1;
    -webkit-transition-property:opacity;
    transition-property:opacity;
    -webkit-transition-duration:200ms;
            transition-duration:200ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-overlay-backdrop.bp3-overlay-exit{
    opacity:1; }
  .bp3-overlay-backdrop.bp3-overlay-exit-active{
    opacity:0;
    -webkit-transition-property:opacity;
    transition-property:opacity;
    -webkit-transition-duration:200ms;
            transition-duration:200ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-overlay-backdrop:focus{
    outline:none; }
  .bp3-overlay-inline .bp3-overlay-backdrop{
    position:absolute; }
.bp3-panel-stack{
  position:relative;
  overflow:hidden; }

.bp3-panel-stack-header{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -ms-flex-negative:0;
      flex-shrink:0;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  z-index:1;
  -webkit-box-shadow:0 1px rgba(16, 22, 26, 0.15);
          box-shadow:0 1px rgba(16, 22, 26, 0.15);
  height:30px; }
  .bp3-dark .bp3-panel-stack-header{
    -webkit-box-shadow:0 1px rgba(255, 255, 255, 0.15);
            box-shadow:0 1px rgba(255, 255, 255, 0.15); }
  .bp3-panel-stack-header > span{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    -webkit-box-flex:1;
        -ms-flex:1;
            flex:1;
    -webkit-box-align:stretch;
        -ms-flex-align:stretch;
            align-items:stretch; }
  .bp3-panel-stack-header .bp3-heading{
    margin:0 5px; }

.bp3-button.bp3-panel-stack-header-back{
  margin-left:5px;
  padding-left:0;
  white-space:nowrap; }
  .bp3-button.bp3-panel-stack-header-back .bp3-icon{
    margin:0 2px; }

.bp3-panel-stack-view{
  position:absolute;
  top:0;
  right:0;
  bottom:0;
  left:0;
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:vertical;
  -webkit-box-direction:normal;
      -ms-flex-direction:column;
          flex-direction:column;
  margin-right:-1px;
  border-right:1px solid rgba(16, 22, 26, 0.15);
  background-color:#ffffff;
  overflow-y:auto; }
  .bp3-dark .bp3-panel-stack-view{
    background-color:#30404d; }

.bp3-panel-stack-push .bp3-panel-stack-enter, .bp3-panel-stack-push .bp3-panel-stack-appear{
  -webkit-transform:translateX(100%);
          transform:translateX(100%);
  opacity:0; }

.bp3-panel-stack-push .bp3-panel-stack-enter-active, .bp3-panel-stack-push .bp3-panel-stack-appear-active{
  -webkit-transform:translate(0%);
          transform:translate(0%);
  opacity:1;
  -webkit-transition-property:opacity, -webkit-transform;
  transition-property:opacity, -webkit-transform;
  transition-property:transform, opacity;
  transition-property:transform, opacity, -webkit-transform;
  -webkit-transition-duration:400ms;
          transition-duration:400ms;
  -webkit-transition-timing-function:ease;
          transition-timing-function:ease;
  -webkit-transition-delay:0;
          transition-delay:0; }

.bp3-panel-stack-push .bp3-panel-stack-exit{
  -webkit-transform:translate(0%);
          transform:translate(0%);
  opacity:1; }

.bp3-panel-stack-push .bp3-panel-stack-exit-active{
  -webkit-transform:translateX(-50%);
          transform:translateX(-50%);
  opacity:0;
  -webkit-transition-property:opacity, -webkit-transform;
  transition-property:opacity, -webkit-transform;
  transition-property:transform, opacity;
  transition-property:transform, opacity, -webkit-transform;
  -webkit-transition-duration:400ms;
          transition-duration:400ms;
  -webkit-transition-timing-function:ease;
          transition-timing-function:ease;
  -webkit-transition-delay:0;
          transition-delay:0; }

.bp3-panel-stack-pop .bp3-panel-stack-enter, .bp3-panel-stack-pop .bp3-panel-stack-appear{
  -webkit-transform:translateX(-50%);
          transform:translateX(-50%);
  opacity:0; }

.bp3-panel-stack-pop .bp3-panel-stack-enter-active, .bp3-panel-stack-pop .bp3-panel-stack-appear-active{
  -webkit-transform:translate(0%);
          transform:translate(0%);
  opacity:1;
  -webkit-transition-property:opacity, -webkit-transform;
  transition-property:opacity, -webkit-transform;
  transition-property:transform, opacity;
  transition-property:transform, opacity, -webkit-transform;
  -webkit-transition-duration:400ms;
          transition-duration:400ms;
  -webkit-transition-timing-function:ease;
          transition-timing-function:ease;
  -webkit-transition-delay:0;
          transition-delay:0; }

.bp3-panel-stack-pop .bp3-panel-stack-exit{
  -webkit-transform:translate(0%);
          transform:translate(0%);
  opacity:1; }

.bp3-panel-stack-pop .bp3-panel-stack-exit-active{
  -webkit-transform:translateX(100%);
          transform:translateX(100%);
  opacity:0;
  -webkit-transition-property:opacity, -webkit-transform;
  transition-property:opacity, -webkit-transform;
  transition-property:transform, opacity;
  transition-property:transform, opacity, -webkit-transform;
  -webkit-transition-duration:400ms;
          transition-duration:400ms;
  -webkit-transition-timing-function:ease;
          transition-timing-function:ease;
  -webkit-transition-delay:0;
          transition-delay:0; }
.bp3-popover{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
  -webkit-transform:scale(1);
          transform:scale(1);
  display:inline-block;
  z-index:20;
  border-radius:3px; }
  .bp3-popover .bp3-popover-arrow{
    position:absolute;
    width:30px;
    height:30px; }
    .bp3-popover .bp3-popover-arrow::before{
      margin:5px;
      width:20px;
      height:20px; }
  .bp3-tether-element-attached-bottom.bp3-tether-target-attached-top > .bp3-popover{
    margin-top:-17px;
    margin-bottom:17px; }
    .bp3-tether-element-attached-bottom.bp3-tether-target-attached-top > .bp3-popover > .bp3-popover-arrow{
      bottom:-11px; }
      .bp3-tether-element-attached-bottom.bp3-tether-target-attached-top > .bp3-popover > .bp3-popover-arrow svg{
        -webkit-transform:rotate(-90deg);
                transform:rotate(-90deg); }
  .bp3-tether-element-attached-left.bp3-tether-target-attached-right > .bp3-popover{
    margin-left:17px; }
    .bp3-tether-element-attached-left.bp3-tether-target-attached-right > .bp3-popover > .bp3-popover-arrow{
      left:-11px; }
      .bp3-tether-element-attached-left.bp3-tether-target-attached-right > .bp3-popover > .bp3-popover-arrow svg{
        -webkit-transform:rotate(0);
                transform:rotate(0); }
  .bp3-tether-element-attached-top.bp3-tether-target-attached-bottom > .bp3-popover{
    margin-top:17px; }
    .bp3-tether-element-attached-top.bp3-tether-target-attached-bottom > .bp3-popover > .bp3-popover-arrow{
      top:-11px; }
      .bp3-tether-element-attached-top.bp3-tether-target-attached-bottom > .bp3-popover > .bp3-popover-arrow svg{
        -webkit-transform:rotate(90deg);
                transform:rotate(90deg); }
  .bp3-tether-element-attached-right.bp3-tether-target-attached-left > .bp3-popover{
    margin-right:17px;
    margin-left:-17px; }
    .bp3-tether-element-attached-right.bp3-tether-target-attached-left > .bp3-popover > .bp3-popover-arrow{
      right:-11px; }
      .bp3-tether-element-attached-right.bp3-tether-target-attached-left > .bp3-popover > .bp3-popover-arrow svg{
        -webkit-transform:rotate(180deg);
                transform:rotate(180deg); }
  .bp3-tether-element-attached-middle > .bp3-popover > .bp3-popover-arrow{
    top:50%;
    -webkit-transform:translateY(-50%);
            transform:translateY(-50%); }
  .bp3-tether-element-attached-center > .bp3-popover > .bp3-popover-arrow{
    right:50%;
    -webkit-transform:translateX(50%);
            transform:translateX(50%); }
  .bp3-tether-element-attached-top.bp3-tether-target-attached-top > .bp3-popover > .bp3-popover-arrow{
    top:-0.3934px; }
  .bp3-tether-element-attached-right.bp3-tether-target-attached-right > .bp3-popover > .bp3-popover-arrow{
    right:-0.3934px; }
  .bp3-tether-element-attached-left.bp3-tether-target-attached-left > .bp3-popover > .bp3-popover-arrow{
    left:-0.3934px; }
  .bp3-tether-element-attached-bottom.bp3-tether-target-attached-bottom > .bp3-popover > .bp3-popover-arrow{
    bottom:-0.3934px; }
  .bp3-tether-element-attached-top.bp3-tether-element-attached-left > .bp3-popover{
    -webkit-transform-origin:top left;
            transform-origin:top left; }
  .bp3-tether-element-attached-top.bp3-tether-element-attached-center > .bp3-popover{
    -webkit-transform-origin:top center;
            transform-origin:top center; }
  .bp3-tether-element-attached-top.bp3-tether-element-attached-right > .bp3-popover{
    -webkit-transform-origin:top right;
            transform-origin:top right; }
  .bp3-tether-element-attached-middle.bp3-tether-element-attached-left > .bp3-popover{
    -webkit-transform-origin:center left;
            transform-origin:center left; }
  .bp3-tether-element-attached-middle.bp3-tether-element-attached-center > .bp3-popover{
    -webkit-transform-origin:center center;
            transform-origin:center center; }
  .bp3-tether-element-attached-middle.bp3-tether-element-attached-right > .bp3-popover{
    -webkit-transform-origin:center right;
            transform-origin:center right; }
  .bp3-tether-element-attached-bottom.bp3-tether-element-attached-left > .bp3-popover{
    -webkit-transform-origin:bottom left;
            transform-origin:bottom left; }
  .bp3-tether-element-attached-bottom.bp3-tether-element-attached-center > .bp3-popover{
    -webkit-transform-origin:bottom center;
            transform-origin:bottom center; }
  .bp3-tether-element-attached-bottom.bp3-tether-element-attached-right > .bp3-popover{
    -webkit-transform-origin:bottom right;
            transform-origin:bottom right; }
  .bp3-popover .bp3-popover-content{
    background:#ffffff;
    color:inherit; }
  .bp3-popover .bp3-popover-arrow::before{
    -webkit-box-shadow:1px 1px 6px rgba(16, 22, 26, 0.2);
            box-shadow:1px 1px 6px rgba(16, 22, 26, 0.2); }
  .bp3-popover .bp3-popover-arrow-border{
    fill:#10161a;
    fill-opacity:0.1; }
  .bp3-popover .bp3-popover-arrow-fill{
    fill:#ffffff; }
  .bp3-popover-enter > .bp3-popover, .bp3-popover-appear > .bp3-popover{
    -webkit-transform:scale(0.3);
            transform:scale(0.3); }
  .bp3-popover-enter-active > .bp3-popover, .bp3-popover-appear-active > .bp3-popover{
    -webkit-transform:scale(1);
            transform:scale(1);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
            transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-popover-exit > .bp3-popover{
    -webkit-transform:scale(1);
            transform:scale(1); }
  .bp3-popover-exit-active > .bp3-popover{
    -webkit-transform:scale(0.3);
            transform:scale(0.3);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
            transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-popover .bp3-popover-content{
    position:relative;
    border-radius:3px; }
  .bp3-popover.bp3-popover-content-sizing .bp3-popover-content{
    max-width:350px;
    padding:20px; }
  .bp3-popover-target + .bp3-overlay .bp3-popover.bp3-popover-content-sizing{
    width:350px; }
  .bp3-popover.bp3-minimal{
    margin:0 !important; }
    .bp3-popover.bp3-minimal .bp3-popover-arrow{
      display:none; }
    .bp3-popover.bp3-minimal.bp3-popover{
      -webkit-transform:scale(1);
              transform:scale(1); }
      .bp3-popover-enter > .bp3-popover.bp3-minimal.bp3-popover, .bp3-popover-appear > .bp3-popover.bp3-minimal.bp3-popover{
        -webkit-transform:scale(1);
                transform:scale(1); }
      .bp3-popover-enter-active > .bp3-popover.bp3-minimal.bp3-popover, .bp3-popover-appear-active > .bp3-popover.bp3-minimal.bp3-popover{
        -webkit-transform:scale(1);
                transform:scale(1);
        -webkit-transition-property:-webkit-transform;
        transition-property:-webkit-transform;
        transition-property:transform;
        transition-property:transform, -webkit-transform;
        -webkit-transition-duration:100ms;
                transition-duration:100ms;
        -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
                transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
        -webkit-transition-delay:0;
                transition-delay:0; }
      .bp3-popover-exit > .bp3-popover.bp3-minimal.bp3-popover{
        -webkit-transform:scale(1);
                transform:scale(1); }
      .bp3-popover-exit-active > .bp3-popover.bp3-minimal.bp3-popover{
        -webkit-transform:scale(1);
                transform:scale(1);
        -webkit-transition-property:-webkit-transform;
        transition-property:-webkit-transform;
        transition-property:transform;
        transition-property:transform, -webkit-transform;
        -webkit-transition-duration:100ms;
                transition-duration:100ms;
        -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
                transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
        -webkit-transition-delay:0;
                transition-delay:0; }
  .bp3-popover.bp3-dark,
  .bp3-dark .bp3-popover{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4); }
    .bp3-popover.bp3-dark .bp3-popover-content,
    .bp3-dark .bp3-popover .bp3-popover-content{
      background:#30404d;
      color:inherit; }
    .bp3-popover.bp3-dark .bp3-popover-arrow::before,
    .bp3-dark .bp3-popover .bp3-popover-arrow::before{
      -webkit-box-shadow:1px 1px 6px rgba(16, 22, 26, 0.4);
              box-shadow:1px 1px 6px rgba(16, 22, 26, 0.4); }
    .bp3-popover.bp3-dark .bp3-popover-arrow-border,
    .bp3-dark .bp3-popover .bp3-popover-arrow-border{
      fill:#10161a;
      fill-opacity:0.2; }
    .bp3-popover.bp3-dark .bp3-popover-arrow-fill,
    .bp3-dark .bp3-popover .bp3-popover-arrow-fill{
      fill:#30404d; }

.bp3-popover-arrow::before{
  display:block;
  position:absolute;
  -webkit-transform:rotate(45deg);
          transform:rotate(45deg);
  border-radius:2px;
  content:""; }

.bp3-tether-pinned .bp3-popover-arrow{
  display:none; }

.bp3-popover-backdrop{
  background:rgba(255, 255, 255, 0); }

.bp3-transition-container{
  opacity:1;
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  z-index:20; }
  .bp3-transition-container.bp3-popover-enter, .bp3-transition-container.bp3-popover-appear{
    opacity:0; }
  .bp3-transition-container.bp3-popover-enter-active, .bp3-transition-container.bp3-popover-appear-active{
    opacity:1;
    -webkit-transition-property:opacity;
    transition-property:opacity;
    -webkit-transition-duration:100ms;
            transition-duration:100ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-transition-container.bp3-popover-exit{
    opacity:1; }
  .bp3-transition-container.bp3-popover-exit-active{
    opacity:0;
    -webkit-transition-property:opacity;
    transition-property:opacity;
    -webkit-transition-duration:100ms;
            transition-duration:100ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-transition-container:focus{
    outline:none; }
  .bp3-transition-container.bp3-popover-leave .bp3-popover-content{
    pointer-events:none; }
  .bp3-transition-container[data-x-out-of-boundaries]{
    display:none; }

span.bp3-popover-target{
  display:inline-block; }

.bp3-popover-wrapper.bp3-fill{
  width:100%; }

.bp3-portal{
  position:absolute;
  top:0;
  right:0;
  left:0; }
@-webkit-keyframes linear-progress-bar-stripes{
  from{
    background-position:0 0; }
  to{
    background-position:30px 0; } }
@keyframes linear-progress-bar-stripes{
  from{
    background-position:0 0; }
  to{
    background-position:30px 0; } }

.bp3-progress-bar{
  display:block;
  position:relative;
  border-radius:40px;
  background:rgba(92, 112, 128, 0.2);
  width:100%;
  height:8px;
  overflow:hidden; }
  .bp3-progress-bar .bp3-progress-meter{
    position:absolute;
    border-radius:40px;
    background:linear-gradient(-45deg, rgba(255, 255, 255, 0.2) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.2) 50%, rgba(255, 255, 255, 0.2) 75%, transparent 75%);
    background-color:rgba(92, 112, 128, 0.8);
    background-size:30px 30px;
    width:100%;
    height:100%;
    -webkit-transition:width 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:width 200ms cubic-bezier(0.4, 1, 0.75, 0.9); }
  .bp3-progress-bar:not(.bp3-no-animation):not(.bp3-no-stripes) .bp3-progress-meter{
    animation:linear-progress-bar-stripes 300ms linear infinite reverse; }
  .bp3-progress-bar.bp3-no-stripes .bp3-progress-meter{
    background-image:none; }

.bp3-dark .bp3-progress-bar{
  background:rgba(16, 22, 26, 0.5); }
  .bp3-dark .bp3-progress-bar .bp3-progress-meter{
    background-color:#8a9ba8; }

.bp3-progress-bar.bp3-intent-primary .bp3-progress-meter{
  background-color:#137cbd; }

.bp3-progress-bar.bp3-intent-success .bp3-progress-meter{
  background-color:#0f9960; }

.bp3-progress-bar.bp3-intent-warning .bp3-progress-meter{
  background-color:#d9822b; }

.bp3-progress-bar.bp3-intent-danger .bp3-progress-meter{
  background-color:#db3737; }
@-webkit-keyframes skeleton-glow{
  from{
    border-color:rgba(206, 217, 224, 0.2);
    background:rgba(206, 217, 224, 0.2); }
  to{
    border-color:rgba(92, 112, 128, 0.2);
    background:rgba(92, 112, 128, 0.2); } }
@keyframes skeleton-glow{
  from{
    border-color:rgba(206, 217, 224, 0.2);
    background:rgba(206, 217, 224, 0.2); }
  to{
    border-color:rgba(92, 112, 128, 0.2);
    background:rgba(92, 112, 128, 0.2); } }
.bp3-skeleton{
  border-color:rgba(206, 217, 224, 0.2) !important;
  border-radius:2px;
  -webkit-box-shadow:none !important;
          box-shadow:none !important;
  background:rgba(206, 217, 224, 0.2);
  background-clip:padding-box !important;
  cursor:default;
  color:transparent !important;
  -webkit-animation:1000ms linear infinite alternate skeleton-glow;
          animation:1000ms linear infinite alternate skeleton-glow;
  pointer-events:none;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-skeleton::before, .bp3-skeleton::after,
  .bp3-skeleton *{
    visibility:hidden !important; }
.bp3-slider{
  width:100%;
  min-width:150px;
  height:40px;
  position:relative;
  outline:none;
  cursor:default;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-slider:hover{
    cursor:pointer; }
  .bp3-slider:active{
    cursor:-webkit-grabbing;
    cursor:grabbing; }
  .bp3-slider.bp3-disabled{
    opacity:0.5;
    cursor:not-allowed; }
  .bp3-slider.bp3-slider-unlabeled{
    height:16px; }

.bp3-slider-track,
.bp3-slider-progress{
  top:5px;
  right:0;
  left:0;
  height:6px;
  position:absolute; }

.bp3-slider-track{
  border-radius:3px;
  overflow:hidden; }

.bp3-slider-progress{
  background:rgba(92, 112, 128, 0.2); }
  .bp3-dark .bp3-slider-progress{
    background:rgba(16, 22, 26, 0.5); }
  .bp3-slider-progress.bp3-intent-primary{
    background-color:#137cbd; }
  .bp3-slider-progress.bp3-intent-success{
    background-color:#0f9960; }
  .bp3-slider-progress.bp3-intent-warning{
    background-color:#d9822b; }
  .bp3-slider-progress.bp3-intent-danger{
    background-color:#db3737; }

.bp3-slider-handle{
  -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
          box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
  background-color:#f5f8fa;
  background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.8)), to(rgba(255, 255, 255, 0)));
  background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.8), rgba(255, 255, 255, 0));
  color:#182026;
  position:absolute;
  top:0;
  left:0;
  border-radius:3px;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.2);
  cursor:pointer;
  width:16px;
  height:16px; }
  .bp3-slider-handle:hover{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-clip:padding-box;
    background-color:#ebf1f5; }
  .bp3-slider-handle:active, .bp3-slider-handle.bp3-active{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#d8e1e8;
    background-image:none; }
  .bp3-slider-handle:disabled, .bp3-slider-handle.bp3-disabled{
    outline:none;
    -webkit-box-shadow:none;
            box-shadow:none;
    background-color:rgba(206, 217, 224, 0.5);
    background-image:none;
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6); }
    .bp3-slider-handle:disabled.bp3-active, .bp3-slider-handle:disabled.bp3-active:hover, .bp3-slider-handle.bp3-disabled.bp3-active, .bp3-slider-handle.bp3-disabled.bp3-active:hover{
      background:rgba(206, 217, 224, 0.7); }
  .bp3-slider-handle:focus{
    z-index:1; }
  .bp3-slider-handle:hover{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 -1px 0 rgba(16, 22, 26, 0.1);
    background-clip:padding-box;
    background-color:#ebf1f5;
    z-index:2;
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 1px 1px rgba(16, 22, 26, 0.2);
    cursor:-webkit-grab;
    cursor:grab; }
  .bp3-slider-handle.bp3-active{
    -webkit-box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
            box-shadow:inset 0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 2px rgba(16, 22, 26, 0.2);
    background-color:#d8e1e8;
    background-image:none;
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 1px rgba(16, 22, 26, 0.1);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), inset 0 1px 1px rgba(16, 22, 26, 0.1);
    cursor:-webkit-grabbing;
    cursor:grabbing; }
  .bp3-disabled .bp3-slider-handle{
    -webkit-box-shadow:none;
            box-shadow:none;
    background:#bfccd6;
    pointer-events:none; }
  .bp3-dark .bp3-slider-handle{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
    background-color:#394b59;
    background-image:-webkit-gradient(linear, left top, left bottom, from(rgba(255, 255, 255, 0.05)), to(rgba(255, 255, 255, 0)));
    background-image:linear-gradient(to bottom, rgba(255, 255, 255, 0.05), rgba(255, 255, 255, 0));
    color:#f5f8fa; }
    .bp3-dark .bp3-slider-handle:hover, .bp3-dark .bp3-slider-handle:active, .bp3-dark .bp3-slider-handle.bp3-active{
      color:#f5f8fa; }
    .bp3-dark .bp3-slider-handle:hover{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.4);
      background-color:#30404d; }
    .bp3-dark .bp3-slider-handle:active, .bp3-dark .bp3-slider-handle.bp3-active{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.6), inset 0 1px 2px rgba(16, 22, 26, 0.2);
      background-color:#202b33;
      background-image:none; }
    .bp3-dark .bp3-slider-handle:disabled, .bp3-dark .bp3-slider-handle.bp3-disabled{
      -webkit-box-shadow:none;
              box-shadow:none;
      background-color:rgba(57, 75, 89, 0.5);
      background-image:none;
      color:rgba(167, 182, 194, 0.6); }
      .bp3-dark .bp3-slider-handle:disabled.bp3-active, .bp3-dark .bp3-slider-handle.bp3-disabled.bp3-active{
        background:rgba(57, 75, 89, 0.7); }
    .bp3-dark .bp3-slider-handle .bp3-button-spinner .bp3-spinner-head{
      background:rgba(16, 22, 26, 0.5);
      stroke:#8a9ba8; }
    .bp3-dark .bp3-slider-handle, .bp3-dark .bp3-slider-handle:hover{
      background-color:#394b59; }
    .bp3-dark .bp3-slider-handle.bp3-active{
      background-color:#293742; }
  .bp3-dark .bp3-disabled .bp3-slider-handle{
    border-color:#5c7080;
    -webkit-box-shadow:none;
            box-shadow:none;
    background:#5c7080; }
  .bp3-slider-handle .bp3-slider-label{
    margin-left:8px;
    border-radius:3px;
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
    background:#394b59;
    color:#f5f8fa; }
    .bp3-dark .bp3-slider-handle .bp3-slider-label{
      -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
      background:#e1e8ed;
      color:#394b59; }
    .bp3-disabled .bp3-slider-handle .bp3-slider-label{
      -webkit-box-shadow:none;
              box-shadow:none; }
  .bp3-slider-handle.bp3-start, .bp3-slider-handle.bp3-end{
    width:8px; }
  .bp3-slider-handle.bp3-start{
    border-top-right-radius:0;
    border-bottom-right-radius:0; }
  .bp3-slider-handle.bp3-end{
    margin-left:8px;
    border-top-left-radius:0;
    border-bottom-left-radius:0; }
    .bp3-slider-handle.bp3-end .bp3-slider-label{
      margin-left:0; }

.bp3-slider-label{
  -webkit-transform:translate(-50%, 20px);
          transform:translate(-50%, 20px);
  display:inline-block;
  position:absolute;
  padding:2px 5px;
  vertical-align:top;
  line-height:1;
  font-size:12px; }

.bp3-slider.bp3-vertical{
  width:40px;
  min-width:40px;
  height:150px; }
  .bp3-slider.bp3-vertical .bp3-slider-track,
  .bp3-slider.bp3-vertical .bp3-slider-progress{
    top:0;
    bottom:0;
    left:5px;
    width:6px;
    height:auto; }
  .bp3-slider.bp3-vertical .bp3-slider-progress{
    top:auto; }
  .bp3-slider.bp3-vertical .bp3-slider-label{
    -webkit-transform:translate(20px, 50%);
            transform:translate(20px, 50%); }
  .bp3-slider.bp3-vertical .bp3-slider-handle{
    top:auto; }
    .bp3-slider.bp3-vertical .bp3-slider-handle .bp3-slider-label{
      margin-top:-8px;
      margin-left:0; }
    .bp3-slider.bp3-vertical .bp3-slider-handle.bp3-end, .bp3-slider.bp3-vertical .bp3-slider-handle.bp3-start{
      margin-left:0;
      width:16px;
      height:8px; }
    .bp3-slider.bp3-vertical .bp3-slider-handle.bp3-start{
      border-top-left-radius:0;
      border-bottom-right-radius:3px; }
      .bp3-slider.bp3-vertical .bp3-slider-handle.bp3-start .bp3-slider-label{
        -webkit-transform:translate(20px);
                transform:translate(20px); }
    .bp3-slider.bp3-vertical .bp3-slider-handle.bp3-end{
      margin-bottom:8px;
      border-top-left-radius:3px;
      border-bottom-left-radius:0;
      border-bottom-right-radius:0; }

@-webkit-keyframes pt-spinner-animation{
  from{
    -webkit-transform:rotate(0deg);
            transform:rotate(0deg); }
  to{
    -webkit-transform:rotate(360deg);
            transform:rotate(360deg); } }

@keyframes pt-spinner-animation{
  from{
    -webkit-transform:rotate(0deg);
            transform:rotate(0deg); }
  to{
    -webkit-transform:rotate(360deg);
            transform:rotate(360deg); } }

.bp3-spinner{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  -webkit-box-pack:center;
      -ms-flex-pack:center;
          justify-content:center;
  overflow:visible;
  vertical-align:middle; }
  .bp3-spinner svg{
    display:block; }
  .bp3-spinner path{
    fill-opacity:0; }
  .bp3-spinner .bp3-spinner-head{
    -webkit-transform-origin:center;
            transform-origin:center;
    -webkit-transition:stroke-dashoffset 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
    transition:stroke-dashoffset 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
    stroke:rgba(92, 112, 128, 0.8);
    stroke-linecap:round; }
  .bp3-spinner .bp3-spinner-track{
    stroke:rgba(92, 112, 128, 0.2); }

.bp3-spinner-animation{
  -webkit-animation:pt-spinner-animation 500ms linear infinite;
          animation:pt-spinner-animation 500ms linear infinite; }
  .bp3-no-spin > .bp3-spinner-animation{
    -webkit-animation:none;
            animation:none; }

.bp3-dark .bp3-spinner .bp3-spinner-head{
  stroke:#8a9ba8; }

.bp3-dark .bp3-spinner .bp3-spinner-track{
  stroke:rgba(16, 22, 26, 0.5); }

.bp3-spinner.bp3-intent-primary .bp3-spinner-head{
  stroke:#137cbd; }

.bp3-spinner.bp3-intent-success .bp3-spinner-head{
  stroke:#0f9960; }

.bp3-spinner.bp3-intent-warning .bp3-spinner-head{
  stroke:#d9822b; }

.bp3-spinner.bp3-intent-danger .bp3-spinner-head{
  stroke:#db3737; }
.bp3-tabs.bp3-vertical{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex; }
  .bp3-tabs.bp3-vertical > .bp3-tab-list{
    -webkit-box-orient:vertical;
    -webkit-box-direction:normal;
        -ms-flex-direction:column;
            flex-direction:column;
    -webkit-box-align:start;
        -ms-flex-align:start;
            align-items:flex-start; }
    .bp3-tabs.bp3-vertical > .bp3-tab-list .bp3-tab{
      border-radius:3px;
      width:100%;
      padding:0 10px; }
      .bp3-tabs.bp3-vertical > .bp3-tab-list .bp3-tab[aria-selected="true"]{
        -webkit-box-shadow:none;
                box-shadow:none;
        background-color:rgba(19, 124, 189, 0.2); }
    .bp3-tabs.bp3-vertical > .bp3-tab-list .bp3-tab-indicator-wrapper .bp3-tab-indicator{
      top:0;
      right:0;
      bottom:0;
      left:0;
      border-radius:3px;
      background-color:rgba(19, 124, 189, 0.2);
      height:auto; }
  .bp3-tabs.bp3-vertical > .bp3-tab-panel{
    margin-top:0;
    padding-left:20px; }

.bp3-tab-list{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  -webkit-box-align:end;
      -ms-flex-align:end;
          align-items:flex-end;
  position:relative;
  margin:0;
  border:none;
  padding:0;
  list-style:none; }
  .bp3-tab-list > *:not(:last-child){
    margin-right:20px; }

.bp3-tab{
  overflow:hidden;
  text-overflow:ellipsis;
  white-space:nowrap;
  word-wrap:normal;
  -webkit-box-flex:0;
      -ms-flex:0 0 auto;
          flex:0 0 auto;
  position:relative;
  cursor:pointer;
  max-width:100%;
  vertical-align:top;
  line-height:30px;
  color:#182026;
  font-size:14px; }
  .bp3-tab a{
    display:block;
    text-decoration:none;
    color:inherit; }
  .bp3-tab-indicator-wrapper ~ .bp3-tab{
    -webkit-box-shadow:none !important;
            box-shadow:none !important;
    background-color:transparent !important; }
  .bp3-tab[aria-disabled="true"]{
    cursor:not-allowed;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-tab[aria-selected="true"]{
    border-radius:0;
    -webkit-box-shadow:inset 0 -3px 0 #106ba3;
            box-shadow:inset 0 -3px 0 #106ba3; }
  .bp3-tab[aria-selected="true"], .bp3-tab:not([aria-disabled="true"]):hover{
    color:#106ba3; }
  .bp3-tab:focus{
    -moz-outline-radius:0; }
  .bp3-large > .bp3-tab{
    line-height:40px;
    font-size:16px; }

.bp3-tab-panel{
  margin-top:20px; }
  .bp3-tab-panel[aria-hidden="true"]{
    display:none; }

.bp3-tab-indicator-wrapper{
  position:absolute;
  top:0;
  left:0;
  -webkit-transform:translateX(0), translateY(0);
          transform:translateX(0), translateY(0);
  -webkit-transition:height, width, -webkit-transform;
  transition:height, width, -webkit-transform;
  transition:height, transform, width;
  transition:height, transform, width, -webkit-transform;
  -webkit-transition-duration:200ms;
          transition-duration:200ms;
  -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
          transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
  pointer-events:none; }
  .bp3-tab-indicator-wrapper .bp3-tab-indicator{
    position:absolute;
    right:0;
    bottom:0;
    left:0;
    background-color:#106ba3;
    height:3px; }
  .bp3-tab-indicator-wrapper.bp3-no-animation{
    -webkit-transition:none;
    transition:none; }

.bp3-dark .bp3-tab{
  color:#f5f8fa; }
  .bp3-dark .bp3-tab[aria-disabled="true"]{
    color:rgba(167, 182, 194, 0.6); }
  .bp3-dark .bp3-tab[aria-selected="true"]{
    -webkit-box-shadow:inset 0 -3px 0 #48aff0;
            box-shadow:inset 0 -3px 0 #48aff0; }
  .bp3-dark .bp3-tab[aria-selected="true"], .bp3-dark .bp3-tab:not([aria-disabled="true"]):hover{
    color:#48aff0; }

.bp3-dark .bp3-tab-indicator{
  background-color:#48aff0; }

.bp3-flex-expander{
  -webkit-box-flex:1;
      -ms-flex:1 1;
          flex:1 1; }
.bp3-tag{
  display:-webkit-inline-box;
  display:-ms-inline-flexbox;
  display:inline-flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  position:relative;
  border:none;
  border-radius:3px;
  -webkit-box-shadow:none;
          box-shadow:none;
  background-color:#5c7080;
  min-width:20px;
  max-width:100%;
  min-height:20px;
  padding:2px 6px;
  line-height:16px;
  color:#f5f8fa;
  font-size:12px; }
  .bp3-tag.bp3-interactive{
    cursor:pointer; }
    .bp3-tag.bp3-interactive:hover{
      background-color:rgba(92, 112, 128, 0.85); }
    .bp3-tag.bp3-interactive.bp3-active, .bp3-tag.bp3-interactive:active{
      background-color:rgba(92, 112, 128, 0.7); }
  .bp3-tag > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-tag > .bp3-fill{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-tag::before,
  .bp3-tag > *{
    margin-right:4px; }
  .bp3-tag:empty::before,
  .bp3-tag > :last-child{
    margin-right:0; }
  .bp3-tag:focus{
    outline:rgba(19, 124, 189, 0.6) auto 2px;
    outline-offset:0;
    -moz-outline-radius:6px; }
  .bp3-tag.bp3-round{
    border-radius:30px;
    padding-right:8px;
    padding-left:8px; }
  .bp3-dark .bp3-tag{
    background-color:#bfccd6;
    color:#182026; }
    .bp3-dark .bp3-tag.bp3-interactive{
      cursor:pointer; }
      .bp3-dark .bp3-tag.bp3-interactive:hover{
        background-color:rgba(191, 204, 214, 0.85); }
      .bp3-dark .bp3-tag.bp3-interactive.bp3-active, .bp3-dark .bp3-tag.bp3-interactive:active{
        background-color:rgba(191, 204, 214, 0.7); }
    .bp3-dark .bp3-tag > .bp3-icon, .bp3-dark .bp3-tag .bp3-icon-standard, .bp3-dark .bp3-tag .bp3-icon-large{
      fill:currentColor; }
  .bp3-tag > .bp3-icon, .bp3-tag .bp3-icon-standard, .bp3-tag .bp3-icon-large{
    fill:#ffffff; }
  .bp3-tag.bp3-large,
  .bp3-large .bp3-tag{
    min-width:30px;
    min-height:30px;
    padding:0 10px;
    line-height:20px;
    font-size:14px; }
    .bp3-tag.bp3-large::before,
    .bp3-tag.bp3-large > *,
    .bp3-large .bp3-tag::before,
    .bp3-large .bp3-tag > *{
      margin-right:7px; }
    .bp3-tag.bp3-large:empty::before,
    .bp3-tag.bp3-large > :last-child,
    .bp3-large .bp3-tag:empty::before,
    .bp3-large .bp3-tag > :last-child{
      margin-right:0; }
    .bp3-tag.bp3-large.bp3-round,
    .bp3-large .bp3-tag.bp3-round{
      padding-right:12px;
      padding-left:12px; }
  .bp3-tag.bp3-intent-primary{
    background:#137cbd;
    color:#ffffff; }
    .bp3-tag.bp3-intent-primary.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-intent-primary.bp3-interactive:hover{
        background-color:rgba(19, 124, 189, 0.85); }
      .bp3-tag.bp3-intent-primary.bp3-interactive.bp3-active, .bp3-tag.bp3-intent-primary.bp3-interactive:active{
        background-color:rgba(19, 124, 189, 0.7); }
  .bp3-tag.bp3-intent-success{
    background:#0f9960;
    color:#ffffff; }
    .bp3-tag.bp3-intent-success.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-intent-success.bp3-interactive:hover{
        background-color:rgba(15, 153, 96, 0.85); }
      .bp3-tag.bp3-intent-success.bp3-interactive.bp3-active, .bp3-tag.bp3-intent-success.bp3-interactive:active{
        background-color:rgba(15, 153, 96, 0.7); }
  .bp3-tag.bp3-intent-warning{
    background:#d9822b;
    color:#ffffff; }
    .bp3-tag.bp3-intent-warning.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-intent-warning.bp3-interactive:hover{
        background-color:rgba(217, 130, 43, 0.85); }
      .bp3-tag.bp3-intent-warning.bp3-interactive.bp3-active, .bp3-tag.bp3-intent-warning.bp3-interactive:active{
        background-color:rgba(217, 130, 43, 0.7); }
  .bp3-tag.bp3-intent-danger{
    background:#db3737;
    color:#ffffff; }
    .bp3-tag.bp3-intent-danger.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-intent-danger.bp3-interactive:hover{
        background-color:rgba(219, 55, 55, 0.85); }
      .bp3-tag.bp3-intent-danger.bp3-interactive.bp3-active, .bp3-tag.bp3-intent-danger.bp3-interactive:active{
        background-color:rgba(219, 55, 55, 0.7); }
  .bp3-tag.bp3-fill{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    width:100%; }
  .bp3-tag.bp3-minimal > .bp3-icon, .bp3-tag.bp3-minimal .bp3-icon-standard, .bp3-tag.bp3-minimal .bp3-icon-large{
    fill:#5c7080; }
  .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]){
    background-color:rgba(138, 155, 168, 0.2);
    color:#182026; }
    .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive:hover{
        background-color:rgba(92, 112, 128, 0.3); }
      .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive.bp3-active, .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive:active{
        background-color:rgba(92, 112, 128, 0.4); }
    .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]){
      color:#f5f8fa; }
      .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive{
        cursor:pointer; }
        .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive:hover{
          background-color:rgba(191, 204, 214, 0.3); }
        .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive.bp3-active, .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]).bp3-interactive:active{
          background-color:rgba(191, 204, 214, 0.4); }
      .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]) > .bp3-icon, .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]) .bp3-icon-standard, .bp3-dark .bp3-tag.bp3-minimal:not([class*="bp3-intent-"]) .bp3-icon-large{
        fill:#a7b6c2; }
  .bp3-tag.bp3-minimal.bp3-intent-primary{
    background-color:rgba(19, 124, 189, 0.15);
    color:#106ba3; }
    .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive:hover{
        background-color:rgba(19, 124, 189, 0.25); }
      .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive.bp3-active, .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive:active{
        background-color:rgba(19, 124, 189, 0.35); }
    .bp3-tag.bp3-minimal.bp3-intent-primary > .bp3-icon, .bp3-tag.bp3-minimal.bp3-intent-primary .bp3-icon-standard, .bp3-tag.bp3-minimal.bp3-intent-primary .bp3-icon-large{
      fill:#137cbd; }
    .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-primary{
      background-color:rgba(19, 124, 189, 0.25);
      color:#48aff0; }
      .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive{
        cursor:pointer; }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive:hover{
          background-color:rgba(19, 124, 189, 0.35); }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive.bp3-active, .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-primary.bp3-interactive:active{
          background-color:rgba(19, 124, 189, 0.45); }
  .bp3-tag.bp3-minimal.bp3-intent-success{
    background-color:rgba(15, 153, 96, 0.15);
    color:#0d8050; }
    .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive:hover{
        background-color:rgba(15, 153, 96, 0.25); }
      .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive.bp3-active, .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive:active{
        background-color:rgba(15, 153, 96, 0.35); }
    .bp3-tag.bp3-minimal.bp3-intent-success > .bp3-icon, .bp3-tag.bp3-minimal.bp3-intent-success .bp3-icon-standard, .bp3-tag.bp3-minimal.bp3-intent-success .bp3-icon-large{
      fill:#0f9960; }
    .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-success{
      background-color:rgba(15, 153, 96, 0.25);
      color:#3dcc91; }
      .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive{
        cursor:pointer; }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive:hover{
          background-color:rgba(15, 153, 96, 0.35); }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive.bp3-active, .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-success.bp3-interactive:active{
          background-color:rgba(15, 153, 96, 0.45); }
  .bp3-tag.bp3-minimal.bp3-intent-warning{
    background-color:rgba(217, 130, 43, 0.15);
    color:#bf7326; }
    .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive:hover{
        background-color:rgba(217, 130, 43, 0.25); }
      .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive.bp3-active, .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive:active{
        background-color:rgba(217, 130, 43, 0.35); }
    .bp3-tag.bp3-minimal.bp3-intent-warning > .bp3-icon, .bp3-tag.bp3-minimal.bp3-intent-warning .bp3-icon-standard, .bp3-tag.bp3-minimal.bp3-intent-warning .bp3-icon-large{
      fill:#d9822b; }
    .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-warning{
      background-color:rgba(217, 130, 43, 0.25);
      color:#ffb366; }
      .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive{
        cursor:pointer; }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive:hover{
          background-color:rgba(217, 130, 43, 0.35); }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive.bp3-active, .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-warning.bp3-interactive:active{
          background-color:rgba(217, 130, 43, 0.45); }
  .bp3-tag.bp3-minimal.bp3-intent-danger{
    background-color:rgba(219, 55, 55, 0.15);
    color:#c23030; }
    .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive{
      cursor:pointer; }
      .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive:hover{
        background-color:rgba(219, 55, 55, 0.25); }
      .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive.bp3-active, .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive:active{
        background-color:rgba(219, 55, 55, 0.35); }
    .bp3-tag.bp3-minimal.bp3-intent-danger > .bp3-icon, .bp3-tag.bp3-minimal.bp3-intent-danger .bp3-icon-standard, .bp3-tag.bp3-minimal.bp3-intent-danger .bp3-icon-large{
      fill:#db3737; }
    .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-danger{
      background-color:rgba(219, 55, 55, 0.25);
      color:#ff7373; }
      .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive{
        cursor:pointer; }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive:hover{
          background-color:rgba(219, 55, 55, 0.35); }
        .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive.bp3-active, .bp3-dark .bp3-tag.bp3-minimal.bp3-intent-danger.bp3-interactive:active{
          background-color:rgba(219, 55, 55, 0.45); }

.bp3-tag-remove{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  opacity:0.5;
  margin-top:-2px;
  margin-right:-6px !important;
  margin-bottom:-2px;
  border:none;
  background:none;
  cursor:pointer;
  padding:2px;
  padding-left:0;
  color:inherit; }
  .bp3-tag-remove:hover{
    opacity:0.8;
    background:none;
    text-decoration:none; }
  .bp3-tag-remove:active{
    opacity:1; }
  .bp3-tag-remove:empty::before{
    line-height:1;
    font-family:"Icons16", sans-serif;
    font-size:16px;
    font-weight:400;
    font-style:normal;
    -moz-osx-font-smoothing:grayscale;
    -webkit-font-smoothing:antialiased;
    content:""; }
  .bp3-large .bp3-tag-remove{
    margin-right:-10px !important;
    padding:5px;
    padding-left:0; }
    .bp3-large .bp3-tag-remove:empty::before{
      line-height:1;
      font-family:"Icons20", sans-serif;
      font-size:20px;
      font-weight:400;
      font-style:normal; }
.bp3-tag-input{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-orient:horizontal;
  -webkit-box-direction:normal;
      -ms-flex-direction:row;
          flex-direction:row;
  -webkit-box-align:start;
      -ms-flex-align:start;
          align-items:flex-start;
  cursor:text;
  height:auto;
  min-height:30px;
  padding-right:0;
  padding-left:5px;
  line-height:inherit; }
  .bp3-tag-input > *{
    -webkit-box-flex:0;
        -ms-flex-positive:0;
            flex-grow:0;
    -ms-flex-negative:0;
        flex-shrink:0; }
  .bp3-tag-input > .bp3-tag-input-values{
    -webkit-box-flex:1;
        -ms-flex-positive:1;
            flex-grow:1;
    -ms-flex-negative:1;
        flex-shrink:1; }
  .bp3-tag-input .bp3-tag-input-icon{
    margin-top:7px;
    margin-right:7px;
    margin-left:2px;
    color:#5c7080; }
  .bp3-tag-input .bp3-tag-input-values{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    -webkit-box-orient:horizontal;
    -webkit-box-direction:normal;
        -ms-flex-direction:row;
            flex-direction:row;
    -ms-flex-wrap:wrap;
        flex-wrap:wrap;
    -webkit-box-align:center;
        -ms-flex-align:center;
            align-items:center;
    -ms-flex-item-align:stretch;
        align-self:stretch;
    margin-top:5px;
    margin-right:7px;
    min-width:0; }
    .bp3-tag-input .bp3-tag-input-values > *{
      -webkit-box-flex:0;
          -ms-flex-positive:0;
              flex-grow:0;
      -ms-flex-negative:0;
          flex-shrink:0; }
    .bp3-tag-input .bp3-tag-input-values > .bp3-fill{
      -webkit-box-flex:1;
          -ms-flex-positive:1;
              flex-grow:1;
      -ms-flex-negative:1;
          flex-shrink:1; }
    .bp3-tag-input .bp3-tag-input-values::before,
    .bp3-tag-input .bp3-tag-input-values > *{
      margin-right:5px; }
    .bp3-tag-input .bp3-tag-input-values:empty::before,
    .bp3-tag-input .bp3-tag-input-values > :last-child{
      margin-right:0; }
    .bp3-tag-input .bp3-tag-input-values:first-child .bp3-input-ghost:first-child{
      padding-left:5px; }
    .bp3-tag-input .bp3-tag-input-values > *{
      margin-bottom:5px; }
  .bp3-tag-input .bp3-tag{
    overflow-wrap:break-word; }
    .bp3-tag-input .bp3-tag.bp3-active{
      outline:rgba(19, 124, 189, 0.6) auto 2px;
      outline-offset:0;
      -moz-outline-radius:6px; }
  .bp3-tag-input .bp3-input-ghost{
    -webkit-box-flex:1;
        -ms-flex:1 1 auto;
            flex:1 1 auto;
    width:80px;
    line-height:20px; }
    .bp3-tag-input .bp3-input-ghost:disabled, .bp3-tag-input .bp3-input-ghost.bp3-disabled{
      cursor:not-allowed; }
  .bp3-tag-input .bp3-button,
  .bp3-tag-input .bp3-spinner{
    margin:3px;
    margin-left:0; }
  .bp3-tag-input .bp3-button{
    min-width:24px;
    min-height:24px;
    padding:0 7px; }
  .bp3-tag-input.bp3-large{
    height:auto;
    min-height:40px; }
    .bp3-tag-input.bp3-large::before,
    .bp3-tag-input.bp3-large > *{
      margin-right:10px; }
    .bp3-tag-input.bp3-large:empty::before,
    .bp3-tag-input.bp3-large > :last-child{
      margin-right:0; }
    .bp3-tag-input.bp3-large .bp3-tag-input-icon{
      margin-top:10px;
      margin-left:5px; }
    .bp3-tag-input.bp3-large .bp3-input-ghost{
      line-height:30px; }
    .bp3-tag-input.bp3-large .bp3-button{
      min-width:30px;
      min-height:30px;
      padding:5px 10px;
      margin:5px;
      margin-left:0; }
    .bp3-tag-input.bp3-large .bp3-spinner{
      margin:8px;
      margin-left:0; }
  .bp3-tag-input.bp3-active{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
    background-color:#ffffff; }
    .bp3-tag-input.bp3-active.bp3-intent-primary{
      -webkit-box-shadow:0 0 0 1px #106ba3, 0 0 0 3px rgba(16, 107, 163, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #106ba3, 0 0 0 3px rgba(16, 107, 163, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-tag-input.bp3-active.bp3-intent-success{
      -webkit-box-shadow:0 0 0 1px #0d8050, 0 0 0 3px rgba(13, 128, 80, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #0d8050, 0 0 0 3px rgba(13, 128, 80, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-tag-input.bp3-active.bp3-intent-warning{
      -webkit-box-shadow:0 0 0 1px #bf7326, 0 0 0 3px rgba(191, 115, 38, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #bf7326, 0 0 0 3px rgba(191, 115, 38, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
    .bp3-tag-input.bp3-active.bp3-intent-danger{
      -webkit-box-shadow:0 0 0 1px #c23030, 0 0 0 3px rgba(194, 48, 48, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2);
              box-shadow:0 0 0 1px #c23030, 0 0 0 3px rgba(194, 48, 48, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.2); }
  .bp3-dark .bp3-tag-input .bp3-tag-input-icon, .bp3-tag-input.bp3-dark .bp3-tag-input-icon{
    color:#a7b6c2; }
  .bp3-dark .bp3-tag-input .bp3-input-ghost, .bp3-tag-input.bp3-dark .bp3-input-ghost{
    color:#f5f8fa; }
    .bp3-dark .bp3-tag-input .bp3-input-ghost::-webkit-input-placeholder, .bp3-tag-input.bp3-dark .bp3-input-ghost::-webkit-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-tag-input .bp3-input-ghost::-moz-placeholder, .bp3-tag-input.bp3-dark .bp3-input-ghost::-moz-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-tag-input .bp3-input-ghost:-ms-input-placeholder, .bp3-tag-input.bp3-dark .bp3-input-ghost:-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-tag-input .bp3-input-ghost::-ms-input-placeholder, .bp3-tag-input.bp3-dark .bp3-input-ghost::-ms-input-placeholder{
      color:rgba(167, 182, 194, 0.6); }
    .bp3-dark .bp3-tag-input .bp3-input-ghost::placeholder, .bp3-tag-input.bp3-dark .bp3-input-ghost::placeholder{
      color:rgba(167, 182, 194, 0.6); }
  .bp3-dark .bp3-tag-input.bp3-active, .bp3-tag-input.bp3-dark.bp3-active{
    -webkit-box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px #137cbd, 0 0 0 1px #137cbd, 0 0 0 3px rgba(19, 124, 189, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
    background-color:rgba(16, 22, 26, 0.3); }
    .bp3-dark .bp3-tag-input.bp3-active.bp3-intent-primary, .bp3-tag-input.bp3-dark.bp3-active.bp3-intent-primary{
      -webkit-box-shadow:0 0 0 1px #106ba3, 0 0 0 3px rgba(16, 107, 163, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #106ba3, 0 0 0 3px rgba(16, 107, 163, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-tag-input.bp3-active.bp3-intent-success, .bp3-tag-input.bp3-dark.bp3-active.bp3-intent-success{
      -webkit-box-shadow:0 0 0 1px #0d8050, 0 0 0 3px rgba(13, 128, 80, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #0d8050, 0 0 0 3px rgba(13, 128, 80, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-tag-input.bp3-active.bp3-intent-warning, .bp3-tag-input.bp3-dark.bp3-active.bp3-intent-warning{
      -webkit-box-shadow:0 0 0 1px #bf7326, 0 0 0 3px rgba(191, 115, 38, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #bf7326, 0 0 0 3px rgba(191, 115, 38, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }
    .bp3-dark .bp3-tag-input.bp3-active.bp3-intent-danger, .bp3-tag-input.bp3-dark.bp3-active.bp3-intent-danger{
      -webkit-box-shadow:0 0 0 1px #c23030, 0 0 0 3px rgba(194, 48, 48, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4);
              box-shadow:0 0 0 1px #c23030, 0 0 0 3px rgba(194, 48, 48, 0.3), inset 0 0 0 1px rgba(16, 22, 26, 0.3), inset 0 1px 1px rgba(16, 22, 26, 0.4); }

.bp3-input-ghost{
  border:none;
  -webkit-box-shadow:none;
          box-shadow:none;
  background:none;
  padding:0; }
  .bp3-input-ghost::-webkit-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input-ghost::-moz-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input-ghost:-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input-ghost::-ms-input-placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input-ghost::placeholder{
    opacity:1;
    color:rgba(92, 112, 128, 0.6); }
  .bp3-input-ghost:focus{
    outline:none !important; }
.bp3-toast{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-align:start;
      -ms-flex-align:start;
          align-items:flex-start;
  position:relative !important;
  margin:20px 0 0;
  border-radius:3px;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
  background-color:#ffffff;
  min-width:300px;
  max-width:500px;
  pointer-events:all; }
  .bp3-toast.bp3-toast-enter, .bp3-toast.bp3-toast-appear{
    -webkit-transform:translateY(-40px);
            transform:translateY(-40px); }
  .bp3-toast.bp3-toast-enter-active, .bp3-toast.bp3-toast-appear-active{
    -webkit-transform:translateY(0);
            transform:translateY(0);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
            transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-toast.bp3-toast-enter ~ .bp3-toast, .bp3-toast.bp3-toast-appear ~ .bp3-toast{
    -webkit-transform:translateY(-40px);
            transform:translateY(-40px); }
  .bp3-toast.bp3-toast-enter-active ~ .bp3-toast, .bp3-toast.bp3-toast-appear-active ~ .bp3-toast{
    -webkit-transform:translateY(0);
            transform:translateY(0);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
            transition-timing-function:cubic-bezier(0.54, 1.12, 0.38, 1.11);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-toast.bp3-toast-exit{
    opacity:1;
    -webkit-filter:blur(0);
            filter:blur(0); }
  .bp3-toast.bp3-toast-exit-active{
    opacity:0;
    -webkit-filter:blur(10px);
            filter:blur(10px);
    -webkit-transition-property:opacity, -webkit-filter;
    transition-property:opacity, -webkit-filter;
    transition-property:opacity, filter;
    transition-property:opacity, filter, -webkit-filter;
    -webkit-transition-duration:300ms;
            transition-duration:300ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-toast.bp3-toast-exit ~ .bp3-toast{
    -webkit-transform:translateY(0);
            transform:translateY(0); }
  .bp3-toast.bp3-toast-exit-active ~ .bp3-toast{
    -webkit-transform:translateY(-40px);
            transform:translateY(-40px);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:100ms;
            transition-duration:100ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:50ms;
            transition-delay:50ms; }
  .bp3-toast .bp3-button-group{
    -webkit-box-flex:0;
        -ms-flex:0 0 auto;
            flex:0 0 auto;
    padding:5px;
    padding-left:0; }
  .bp3-toast > .bp3-icon{
    margin:12px;
    margin-right:0;
    color:#5c7080; }
  .bp3-toast.bp3-dark,
  .bp3-dark .bp3-toast{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
    background-color:#394b59; }
    .bp3-toast.bp3-dark > .bp3-icon,
    .bp3-dark .bp3-toast > .bp3-icon{
      color:#a7b6c2; }
  .bp3-toast[class*="bp3-intent-"] a{
    color:rgba(255, 255, 255, 0.7); }
    .bp3-toast[class*="bp3-intent-"] a:hover{
      color:#ffffff; }
  .bp3-toast[class*="bp3-intent-"] > .bp3-icon{
    color:#ffffff; }
  .bp3-toast[class*="bp3-intent-"] .bp3-button, .bp3-toast[class*="bp3-intent-"] .bp3-button::before,
  .bp3-toast[class*="bp3-intent-"] .bp3-button .bp3-icon, .bp3-toast[class*="bp3-intent-"] .bp3-button:active{
    color:rgba(255, 255, 255, 0.7) !important; }
  .bp3-toast[class*="bp3-intent-"] .bp3-button:focus{
    outline-color:rgba(255, 255, 255, 0.5); }
  .bp3-toast[class*="bp3-intent-"] .bp3-button:hover{
    background-color:rgba(255, 255, 255, 0.15) !important;
    color:#ffffff !important; }
  .bp3-toast[class*="bp3-intent-"] .bp3-button:active{
    background-color:rgba(255, 255, 255, 0.3) !important;
    color:#ffffff !important; }
  .bp3-toast[class*="bp3-intent-"] .bp3-button::after{
    background:rgba(255, 255, 255, 0.3) !important; }
  .bp3-toast.bp3-intent-primary{
    background-color:#137cbd;
    color:#ffffff; }
  .bp3-toast.bp3-intent-success{
    background-color:#0f9960;
    color:#ffffff; }
  .bp3-toast.bp3-intent-warning{
    background-color:#d9822b;
    color:#ffffff; }
  .bp3-toast.bp3-intent-danger{
    background-color:#db3737;
    color:#ffffff; }

.bp3-toast-message{
  -webkit-box-flex:1;
      -ms-flex:1 1 auto;
          flex:1 1 auto;
  padding:11px;
  word-break:break-word; }

.bp3-toast-container{
  display:-webkit-box !important;
  display:-ms-flexbox !important;
  display:flex !important;
  -webkit-box-orient:vertical;
  -webkit-box-direction:normal;
      -ms-flex-direction:column;
          flex-direction:column;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  position:fixed;
  right:0;
  left:0;
  z-index:40;
  overflow:hidden;
  padding:0 20px 20px;
  pointer-events:none; }
  .bp3-toast-container.bp3-toast-container-top{
    top:0;
    bottom:auto; }
  .bp3-toast-container.bp3-toast-container-bottom{
    -webkit-box-orient:vertical;
    -webkit-box-direction:reverse;
        -ms-flex-direction:column-reverse;
            flex-direction:column-reverse;
    top:auto;
    bottom:0; }
  .bp3-toast-container.bp3-toast-container-left{
    -webkit-box-align:start;
        -ms-flex-align:start;
            align-items:flex-start; }
  .bp3-toast-container.bp3-toast-container-right{
    -webkit-box-align:end;
        -ms-flex-align:end;
            align-items:flex-end; }

.bp3-toast-container-bottom .bp3-toast.bp3-toast-enter:not(.bp3-toast-enter-active),
.bp3-toast-container-bottom .bp3-toast.bp3-toast-enter:not(.bp3-toast-enter-active) ~ .bp3-toast, .bp3-toast-container-bottom .bp3-toast.bp3-toast-appear:not(.bp3-toast-appear-active),
.bp3-toast-container-bottom .bp3-toast.bp3-toast-appear:not(.bp3-toast-appear-active) ~ .bp3-toast,
.bp3-toast-container-bottom .bp3-toast.bp3-toast-leave-active ~ .bp3-toast{
  -webkit-transform:translateY(60px);
          transform:translateY(60px); }
.bp3-tooltip{
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 2px 4px rgba(16, 22, 26, 0.2), 0 8px 24px rgba(16, 22, 26, 0.2);
  -webkit-transform:scale(1);
          transform:scale(1); }
  .bp3-tooltip .bp3-popover-arrow{
    position:absolute;
    width:22px;
    height:22px; }
    .bp3-tooltip .bp3-popover-arrow::before{
      margin:4px;
      width:14px;
      height:14px; }
  .bp3-tether-element-attached-bottom.bp3-tether-target-attached-top > .bp3-tooltip{
    margin-top:-11px;
    margin-bottom:11px; }
    .bp3-tether-element-attached-bottom.bp3-tether-target-attached-top > .bp3-tooltip > .bp3-popover-arrow{
      bottom:-8px; }
      .bp3-tether-element-attached-bottom.bp3-tether-target-attached-top > .bp3-tooltip > .bp3-popover-arrow svg{
        -webkit-transform:rotate(-90deg);
                transform:rotate(-90deg); }
  .bp3-tether-element-attached-left.bp3-tether-target-attached-right > .bp3-tooltip{
    margin-left:11px; }
    .bp3-tether-element-attached-left.bp3-tether-target-attached-right > .bp3-tooltip > .bp3-popover-arrow{
      left:-8px; }
      .bp3-tether-element-attached-left.bp3-tether-target-attached-right > .bp3-tooltip > .bp3-popover-arrow svg{
        -webkit-transform:rotate(0);
                transform:rotate(0); }
  .bp3-tether-element-attached-top.bp3-tether-target-attached-bottom > .bp3-tooltip{
    margin-top:11px; }
    .bp3-tether-element-attached-top.bp3-tether-target-attached-bottom > .bp3-tooltip > .bp3-popover-arrow{
      top:-8px; }
      .bp3-tether-element-attached-top.bp3-tether-target-attached-bottom > .bp3-tooltip > .bp3-popover-arrow svg{
        -webkit-transform:rotate(90deg);
                transform:rotate(90deg); }
  .bp3-tether-element-attached-right.bp3-tether-target-attached-left > .bp3-tooltip{
    margin-right:11px;
    margin-left:-11px; }
    .bp3-tether-element-attached-right.bp3-tether-target-attached-left > .bp3-tooltip > .bp3-popover-arrow{
      right:-8px; }
      .bp3-tether-element-attached-right.bp3-tether-target-attached-left > .bp3-tooltip > .bp3-popover-arrow svg{
        -webkit-transform:rotate(180deg);
                transform:rotate(180deg); }
  .bp3-tether-element-attached-middle > .bp3-tooltip > .bp3-popover-arrow{
    top:50%;
    -webkit-transform:translateY(-50%);
            transform:translateY(-50%); }
  .bp3-tether-element-attached-center > .bp3-tooltip > .bp3-popover-arrow{
    right:50%;
    -webkit-transform:translateX(50%);
            transform:translateX(50%); }
  .bp3-tether-element-attached-top.bp3-tether-target-attached-top > .bp3-tooltip > .bp3-popover-arrow{
    top:-0.22183px; }
  .bp3-tether-element-attached-right.bp3-tether-target-attached-right > .bp3-tooltip > .bp3-popover-arrow{
    right:-0.22183px; }
  .bp3-tether-element-attached-left.bp3-tether-target-attached-left > .bp3-tooltip > .bp3-popover-arrow{
    left:-0.22183px; }
  .bp3-tether-element-attached-bottom.bp3-tether-target-attached-bottom > .bp3-tooltip > .bp3-popover-arrow{
    bottom:-0.22183px; }
  .bp3-tether-element-attached-top.bp3-tether-element-attached-left > .bp3-tooltip{
    -webkit-transform-origin:top left;
            transform-origin:top left; }
  .bp3-tether-element-attached-top.bp3-tether-element-attached-center > .bp3-tooltip{
    -webkit-transform-origin:top center;
            transform-origin:top center; }
  .bp3-tether-element-attached-top.bp3-tether-element-attached-right > .bp3-tooltip{
    -webkit-transform-origin:top right;
            transform-origin:top right; }
  .bp3-tether-element-attached-middle.bp3-tether-element-attached-left > .bp3-tooltip{
    -webkit-transform-origin:center left;
            transform-origin:center left; }
  .bp3-tether-element-attached-middle.bp3-tether-element-attached-center > .bp3-tooltip{
    -webkit-transform-origin:center center;
            transform-origin:center center; }
  .bp3-tether-element-attached-middle.bp3-tether-element-attached-right > .bp3-tooltip{
    -webkit-transform-origin:center right;
            transform-origin:center right; }
  .bp3-tether-element-attached-bottom.bp3-tether-element-attached-left > .bp3-tooltip{
    -webkit-transform-origin:bottom left;
            transform-origin:bottom left; }
  .bp3-tether-element-attached-bottom.bp3-tether-element-attached-center > .bp3-tooltip{
    -webkit-transform-origin:bottom center;
            transform-origin:bottom center; }
  .bp3-tether-element-attached-bottom.bp3-tether-element-attached-right > .bp3-tooltip{
    -webkit-transform-origin:bottom right;
            transform-origin:bottom right; }
  .bp3-tooltip .bp3-popover-content{
    background:#394b59;
    color:#f5f8fa; }
  .bp3-tooltip .bp3-popover-arrow::before{
    -webkit-box-shadow:1px 1px 6px rgba(16, 22, 26, 0.2);
            box-shadow:1px 1px 6px rgba(16, 22, 26, 0.2); }
  .bp3-tooltip .bp3-popover-arrow-border{
    fill:#10161a;
    fill-opacity:0.1; }
  .bp3-tooltip .bp3-popover-arrow-fill{
    fill:#394b59; }
  .bp3-popover-enter > .bp3-tooltip, .bp3-popover-appear > .bp3-tooltip{
    -webkit-transform:scale(0.8);
            transform:scale(0.8); }
  .bp3-popover-enter-active > .bp3-tooltip, .bp3-popover-appear-active > .bp3-tooltip{
    -webkit-transform:scale(1);
            transform:scale(1);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:100ms;
            transition-duration:100ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-popover-exit > .bp3-tooltip{
    -webkit-transform:scale(1);
            transform:scale(1); }
  .bp3-popover-exit-active > .bp3-tooltip{
    -webkit-transform:scale(0.8);
            transform:scale(0.8);
    -webkit-transition-property:-webkit-transform;
    transition-property:-webkit-transform;
    transition-property:transform;
    transition-property:transform, -webkit-transform;
    -webkit-transition-duration:100ms;
            transition-duration:100ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-tooltip .bp3-popover-content{
    padding:10px 12px; }
  .bp3-tooltip.bp3-dark,
  .bp3-dark .bp3-tooltip{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 2px 4px rgba(16, 22, 26, 0.4), 0 8px 24px rgba(16, 22, 26, 0.4); }
    .bp3-tooltip.bp3-dark .bp3-popover-content,
    .bp3-dark .bp3-tooltip .bp3-popover-content{
      background:#e1e8ed;
      color:#394b59; }
    .bp3-tooltip.bp3-dark .bp3-popover-arrow::before,
    .bp3-dark .bp3-tooltip .bp3-popover-arrow::before{
      -webkit-box-shadow:1px 1px 6px rgba(16, 22, 26, 0.4);
              box-shadow:1px 1px 6px rgba(16, 22, 26, 0.4); }
    .bp3-tooltip.bp3-dark .bp3-popover-arrow-border,
    .bp3-dark .bp3-tooltip .bp3-popover-arrow-border{
      fill:#10161a;
      fill-opacity:0.2; }
    .bp3-tooltip.bp3-dark .bp3-popover-arrow-fill,
    .bp3-dark .bp3-tooltip .bp3-popover-arrow-fill{
      fill:#e1e8ed; }
  .bp3-tooltip.bp3-intent-primary .bp3-popover-content{
    background:#137cbd;
    color:#ffffff; }
  .bp3-tooltip.bp3-intent-primary .bp3-popover-arrow-fill{
    fill:#137cbd; }
  .bp3-tooltip.bp3-intent-success .bp3-popover-content{
    background:#0f9960;
    color:#ffffff; }
  .bp3-tooltip.bp3-intent-success .bp3-popover-arrow-fill{
    fill:#0f9960; }
  .bp3-tooltip.bp3-intent-warning .bp3-popover-content{
    background:#d9822b;
    color:#ffffff; }
  .bp3-tooltip.bp3-intent-warning .bp3-popover-arrow-fill{
    fill:#d9822b; }
  .bp3-tooltip.bp3-intent-danger .bp3-popover-content{
    background:#db3737;
    color:#ffffff; }
  .bp3-tooltip.bp3-intent-danger .bp3-popover-arrow-fill{
    fill:#db3737; }

.bp3-tooltip-indicator{
  border-bottom:dotted 1px;
  cursor:help; }
.bp3-tree .bp3-icon, .bp3-tree .bp3-icon-standard, .bp3-tree .bp3-icon-large{
  color:#5c7080; }
  .bp3-tree .bp3-icon.bp3-intent-primary, .bp3-tree .bp3-icon-standard.bp3-intent-primary, .bp3-tree .bp3-icon-large.bp3-intent-primary{
    color:#137cbd; }
  .bp3-tree .bp3-icon.bp3-intent-success, .bp3-tree .bp3-icon-standard.bp3-intent-success, .bp3-tree .bp3-icon-large.bp3-intent-success{
    color:#0f9960; }
  .bp3-tree .bp3-icon.bp3-intent-warning, .bp3-tree .bp3-icon-standard.bp3-intent-warning, .bp3-tree .bp3-icon-large.bp3-intent-warning{
    color:#d9822b; }
  .bp3-tree .bp3-icon.bp3-intent-danger, .bp3-tree .bp3-icon-standard.bp3-intent-danger, .bp3-tree .bp3-icon-large.bp3-intent-danger{
    color:#db3737; }

.bp3-tree-node-list{
  margin:0;
  padding-left:0;
  list-style:none; }

.bp3-tree-root{
  position:relative;
  background-color:transparent;
  cursor:default;
  padding-left:0; }

.bp3-tree-node-content-0{
  padding-left:0px; }

.bp3-tree-node-content-1{
  padding-left:23px; }

.bp3-tree-node-content-2{
  padding-left:46px; }

.bp3-tree-node-content-3{
  padding-left:69px; }

.bp3-tree-node-content-4{
  padding-left:92px; }

.bp3-tree-node-content-5{
  padding-left:115px; }

.bp3-tree-node-content-6{
  padding-left:138px; }

.bp3-tree-node-content-7{
  padding-left:161px; }

.bp3-tree-node-content-8{
  padding-left:184px; }

.bp3-tree-node-content-9{
  padding-left:207px; }

.bp3-tree-node-content-10{
  padding-left:230px; }

.bp3-tree-node-content-11{
  padding-left:253px; }

.bp3-tree-node-content-12{
  padding-left:276px; }

.bp3-tree-node-content-13{
  padding-left:299px; }

.bp3-tree-node-content-14{
  padding-left:322px; }

.bp3-tree-node-content-15{
  padding-left:345px; }

.bp3-tree-node-content-16{
  padding-left:368px; }

.bp3-tree-node-content-17{
  padding-left:391px; }

.bp3-tree-node-content-18{
  padding-left:414px; }

.bp3-tree-node-content-19{
  padding-left:437px; }

.bp3-tree-node-content-20{
  padding-left:460px; }

.bp3-tree-node-content{
  display:-webkit-box;
  display:-ms-flexbox;
  display:flex;
  -webkit-box-align:center;
      -ms-flex-align:center;
          align-items:center;
  width:100%;
  height:30px;
  padding-right:5px; }
  .bp3-tree-node-content:hover{
    background-color:rgba(191, 204, 214, 0.4); }

.bp3-tree-node-caret,
.bp3-tree-node-caret-none{
  min-width:30px; }

.bp3-tree-node-caret{
  color:#5c7080;
  -webkit-transform:rotate(0deg);
          transform:rotate(0deg);
  cursor:pointer;
  padding:7px;
  -webkit-transition:-webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:-webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9);
  transition:transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9), -webkit-transform 200ms cubic-bezier(0.4, 1, 0.75, 0.9); }
  .bp3-tree-node-caret:hover{
    color:#182026; }
  .bp3-dark .bp3-tree-node-caret{
    color:#a7b6c2; }
    .bp3-dark .bp3-tree-node-caret:hover{
      color:#f5f8fa; }
  .bp3-tree-node-caret.bp3-tree-node-caret-open{
    -webkit-transform:rotate(90deg);
            transform:rotate(90deg); }
  .bp3-tree-node-caret.bp3-icon-standard::before{
    content:""; }

.bp3-tree-node-icon{
  position:relative;
  margin-right:7px; }

.bp3-tree-node-label{
  overflow:hidden;
  text-overflow:ellipsis;
  white-space:nowrap;
  word-wrap:normal;
  -webkit-box-flex:1;
      -ms-flex:1 1 auto;
          flex:1 1 auto;
  position:relative;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-tree-node-label span{
    display:inline; }

.bp3-tree-node-secondary-label{
  padding:0 5px;
  -webkit-user-select:none;
     -moz-user-select:none;
      -ms-user-select:none;
          user-select:none; }
  .bp3-tree-node-secondary-label .bp3-popover-wrapper,
  .bp3-tree-node-secondary-label .bp3-popover-target{
    display:-webkit-box;
    display:-ms-flexbox;
    display:flex;
    -webkit-box-align:center;
        -ms-flex-align:center;
            align-items:center; }

.bp3-tree-node.bp3-disabled .bp3-tree-node-content{
  background-color:inherit;
  cursor:not-allowed;
  color:rgba(92, 112, 128, 0.6); }

.bp3-tree-node.bp3-disabled .bp3-tree-node-caret,
.bp3-tree-node.bp3-disabled .bp3-tree-node-icon{
  cursor:not-allowed;
  color:rgba(92, 112, 128, 0.6); }

.bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content{
  background-color:#137cbd; }
  .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content,
  .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content .bp3-icon, .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content .bp3-icon-standard, .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content .bp3-icon-large{
    color:#ffffff; }
  .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content .bp3-tree-node-caret::before{
    color:rgba(255, 255, 255, 0.7); }
  .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content .bp3-tree-node-caret:hover::before{
    color:#ffffff; }

.bp3-dark .bp3-tree-node-content:hover{
  background-color:rgba(92, 112, 128, 0.3); }

.bp3-dark .bp3-tree .bp3-icon, .bp3-dark .bp3-tree .bp3-icon-standard, .bp3-dark .bp3-tree .bp3-icon-large{
  color:#a7b6c2; }
  .bp3-dark .bp3-tree .bp3-icon.bp3-intent-primary, .bp3-dark .bp3-tree .bp3-icon-standard.bp3-intent-primary, .bp3-dark .bp3-tree .bp3-icon-large.bp3-intent-primary{
    color:#137cbd; }
  .bp3-dark .bp3-tree .bp3-icon.bp3-intent-success, .bp3-dark .bp3-tree .bp3-icon-standard.bp3-intent-success, .bp3-dark .bp3-tree .bp3-icon-large.bp3-intent-success{
    color:#0f9960; }
  .bp3-dark .bp3-tree .bp3-icon.bp3-intent-warning, .bp3-dark .bp3-tree .bp3-icon-standard.bp3-intent-warning, .bp3-dark .bp3-tree .bp3-icon-large.bp3-intent-warning{
    color:#d9822b; }
  .bp3-dark .bp3-tree .bp3-icon.bp3-intent-danger, .bp3-dark .bp3-tree .bp3-icon-standard.bp3-intent-danger, .bp3-dark .bp3-tree .bp3-icon-large.bp3-intent-danger{
    color:#db3737; }

.bp3-dark .bp3-tree-node.bp3-tree-node-selected > .bp3-tree-node-content{
  background-color:#137cbd; }
/*!

Copyright 2017-present Palantir Technologies, Inc. All rights reserved.
Licensed under the Apache License, Version 2.0.

*/
.bp3-omnibar{
  -webkit-filter:blur(0);
          filter:blur(0);
  opacity:1;
  top:20vh;
  left:calc(50% - 250px);
  z-index:21;
  border-radius:3px;
  -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
          box-shadow:0 0 0 1px rgba(16, 22, 26, 0.1), 0 4px 8px rgba(16, 22, 26, 0.2), 0 18px 46px 6px rgba(16, 22, 26, 0.2);
  background-color:#ffffff;
  width:500px; }
  .bp3-omnibar.bp3-overlay-enter, .bp3-omnibar.bp3-overlay-appear{
    -webkit-filter:blur(20px);
            filter:blur(20px);
    opacity:0.2; }
  .bp3-omnibar.bp3-overlay-enter-active, .bp3-omnibar.bp3-overlay-appear-active{
    -webkit-filter:blur(0);
            filter:blur(0);
    opacity:1;
    -webkit-transition-property:opacity, -webkit-filter;
    transition-property:opacity, -webkit-filter;
    transition-property:filter, opacity;
    transition-property:filter, opacity, -webkit-filter;
    -webkit-transition-duration:200ms;
            transition-duration:200ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-omnibar.bp3-overlay-exit{
    -webkit-filter:blur(0);
            filter:blur(0);
    opacity:1; }
  .bp3-omnibar.bp3-overlay-exit-active{
    -webkit-filter:blur(20px);
            filter:blur(20px);
    opacity:0.2;
    -webkit-transition-property:opacity, -webkit-filter;
    transition-property:opacity, -webkit-filter;
    transition-property:filter, opacity;
    transition-property:filter, opacity, -webkit-filter;
    -webkit-transition-duration:200ms;
            transition-duration:200ms;
    -webkit-transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
            transition-timing-function:cubic-bezier(0.4, 1, 0.75, 0.9);
    -webkit-transition-delay:0;
            transition-delay:0; }
  .bp3-omnibar .bp3-input{
    border-radius:0;
    background-color:transparent; }
    .bp3-omnibar .bp3-input, .bp3-omnibar .bp3-input:focus{
      -webkit-box-shadow:none;
              box-shadow:none; }
  .bp3-omnibar .bp3-menu{
    border-radius:0;
    -webkit-box-shadow:inset 0 1px 0 rgba(16, 22, 26, 0.15);
            box-shadow:inset 0 1px 0 rgba(16, 22, 26, 0.15);
    background-color:transparent;
    max-height:calc(60vh - 40px);
    overflow:auto; }
    .bp3-omnibar .bp3-menu:empty{
      display:none; }
  .bp3-dark .bp3-omnibar, .bp3-omnibar.bp3-dark{
    -webkit-box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
            box-shadow:0 0 0 1px rgba(16, 22, 26, 0.2), 0 4px 8px rgba(16, 22, 26, 0.4), 0 18px 46px 6px rgba(16, 22, 26, 0.4);
    background-color:#30404d; }

.bp3-omnibar-overlay .bp3-overlay-backdrop{
  background-color:rgba(16, 22, 26, 0.2); }

.bp3-select-popover .bp3-popover-content{
  padding:5px; }

.bp3-select-popover .bp3-input-group{
  margin-bottom:0; }

.bp3-select-popover .bp3-menu{
  max-width:400px;
  max-height:300px;
  overflow:auto;
  padding:0; }
  .bp3-select-popover .bp3-menu:not(:first-child){
    padding-top:5px; }

.bp3-multi-select{
  min-width:150px; }

.bp3-multi-select-popover .bp3-menu{
  max-width:400px;
  max-height:300px;
  overflow:auto; }

.bp3-select-popover .bp3-popover-content{
  padding:5px; }

.bp3-select-popover .bp3-input-group{
  margin-bottom:0; }

.bp3-select-popover .bp3-menu{
  max-width:400px;
  max-height:300px;
  overflow:auto;
  padding:0; }
  .bp3-select-popover .bp3-menu:not(:first-child){
    padding-top:5px; }
/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensureUiComponents() in @jupyterlab/buildutils */

/**
 * (DEPRECATED) Support for consuming icons as CSS background images
 */

/* Icons urls */

:root {
  --jp-icon-add: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTE5IDEzaC02djZoLTJ2LTZINXYtMmg2VjVoMnY2aDZ2MnoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-bug: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTIwIDhoLTIuODFjLS40NS0uNzgtMS4wNy0xLjQ1LTEuODItMS45NkwxNyA0LjQxIDE1LjU5IDNsLTIuMTcgMi4xN0MxMi45NiA1LjA2IDEyLjQ5IDUgMTIgNWMtLjQ5IDAtLjk2LjA2LTEuNDEuMTdMOC40MSAzIDcgNC40MWwxLjYyIDEuNjNDNy44OCA2LjU1IDcuMjYgNy4yMiA2LjgxIDhINHYyaDIuMDljLS4wNS4zMy0uMDkuNjYtLjA5IDF2MUg0djJoMnYxYzAgLjM0LjA0LjY3LjA5IDFINHYyaDIuODFjMS4wNCAxLjc5IDIuOTcgMyA1LjE5IDNzNC4xNS0xLjIxIDUuMTktM0gyMHYtMmgtMi4wOWMuMDUtLjMzLjA5LS42Ni4wOS0xdi0xaDJ2LTJoLTJ2LTFjMC0uMzQtLjA0LS42Ny0uMDktMUgyMFY4em0tNiA4aC00di0yaDR2MnptMC00aC00di0yaDR2MnoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-build: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTYiIHZpZXdCb3g9IjAgMCAyNCAyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTE0LjkgMTcuNDVDMTYuMjUgMTcuNDUgMTcuMzUgMTYuMzUgMTcuMzUgMTVDMTcuMzUgMTMuNjUgMTYuMjUgMTIuNTUgMTQuOSAxMi41NUMxMy41NCAxMi41NSAxMi40NSAxMy42NSAxMi40NSAxNUMxMi40NSAxNi4zNSAxMy41NCAxNy40NSAxNC45IDE3LjQ1Wk0yMC4xIDE1LjY4TDIxLjU4IDE2Ljg0QzIxLjcxIDE2Ljk1IDIxLjc1IDE3LjEzIDIxLjY2IDE3LjI5TDIwLjI2IDE5LjcxQzIwLjE3IDE5Ljg2IDIwIDE5LjkyIDE5LjgzIDE5Ljg2TDE4LjA5IDE5LjE2QzE3LjczIDE5LjQ0IDE3LjMzIDE5LjY3IDE2LjkxIDE5Ljg1TDE2LjY0IDIxLjdDMTYuNjIgMjEuODcgMTYuNDcgMjIgMTYuMyAyMkgxMy41QzEzLjMyIDIyIDEzLjE4IDIxLjg3IDEzLjE1IDIxLjdMMTIuODkgMTkuODVDMTIuNDYgMTkuNjcgMTIuMDcgMTkuNDQgMTEuNzEgMTkuMTZMOS45NjAwMiAxOS44NkM5LjgxMDAyIDE5LjkyIDkuNjIwMDIgMTkuODYgOS41NDAwMiAxOS43MUw4LjE0MDAyIDE3LjI5QzguMDUwMDIgMTcuMTMgOC4wOTAwMiAxNi45NSA4LjIyMDAyIDE2Ljg0TDkuNzAwMDIgMTUuNjhMOS42NTAwMSAxNUw5LjcwMDAyIDE0LjMxTDguMjIwMDIgMTMuMTZDOC4wOTAwMiAxMy4wNSA4LjA1MDAyIDEyLjg2IDguMTQwMDIgMTIuNzFMOS41NDAwMiAxMC4yOUM5LjYyMDAyIDEwLjEzIDkuODEwMDIgMTAuMDcgOS45NjAwMiAxMC4xM0wxMS43MSAxMC44NEMxMi4wNyAxMC41NiAxMi40NiAxMC4zMiAxMi44OSAxMC4xNUwxMy4xNSA4LjI4OTk4QzEzLjE4IDguMTI5OTggMTMuMzIgNy45OTk5OCAxMy41IDcuOTk5OThIMTYuM0MxNi40NyA3Ljk5OTk4IDE2LjYyIDguMTI5OTggMTYuNjQgOC4yODk5OEwxNi45MSAxMC4xNUMxNy4zMyAxMC4zMiAxNy43MyAxMC41NiAxOC4wOSAxMC44NEwxOS44MyAxMC4xM0MyMCAxMC4wNyAyMC4xNyAxMC4xMyAyMC4yNiAxMC4yOUwyMS42NiAxMi43MUMyMS43NSAxMi44NiAyMS43MSAxMy4wNSAyMS41OCAxMy4xNkwyMC4xIDE0LjMxTDIwLjE1IDE1TDIwLjEgMTUuNjhaIi8+CiAgICA8cGF0aCBkPSJNNy4zMjk2NiA3LjQ0NDU0QzguMDgzMSA3LjAwOTU0IDguMzM5MzIgNi4wNTMzMiA3LjkwNDMyIDUuMjk5ODhDNy40NjkzMiA0LjU0NjQzIDYuNTA4MSA0LjI4MTU2IDUuNzU0NjYgNC43MTY1NkM1LjM5MTc2IDQuOTI2MDggNS4xMjY5NSA1LjI3MTE4IDUuMDE4NDkgNS42NzU5NEM0LjkxMDA0IDYuMDgwNzEgNC45NjY4MiA2LjUxMTk4IDUuMTc2MzQgNi44NzQ4OEM1LjYxMTM0IDcuNjI4MzIgNi41NzYyMiA3Ljg3OTU0IDcuMzI5NjYgNy40NDQ1NFpNOS42NTcxOCA0Ljc5NTkzTDEwLjg2NzIgNC45NTE3OUMxMC45NjI4IDQuOTc3NDEgMTEuMDQwMiA1LjA3MTMzIDExLjAzODIgNS4xODc5M0wxMS4wMzg4IDYuOTg4OTNDMTEuMDQ1NSA3LjEwMDU0IDEwLjk2MTYgNy4xOTUxOCAxMC44NTUgNy4yMTA1NEw5LjY2MDAxIDcuMzgwODNMOS4yMzkxNSA4LjEzMTg4TDkuNjY5NjEgOS4yNTc0NUM5LjcwNzI5IDkuMzYyNzEgOS42NjkzNCA5LjQ3Njk5IDkuNTc0MDggOS41MzE5OUw4LjAxNTIzIDEwLjQzMkM3LjkxMTMxIDEwLjQ5MiA3Ljc5MzM3IDEwLjQ2NzcgNy43MjEwNSAxMC4zODI0TDYuOTg3NDggOS40MzE4OEw2LjEwOTMxIDkuNDMwODNMNS4zNDcwNCAxMC4zOTA1QzUuMjg5MDkgMTAuNDcwMiA1LjE3MzgzIDEwLjQ5MDUgNS4wNzE4NyAxMC40MzM5TDMuNTEyNDUgOS41MzI5M0MzLjQxMDQ5IDkuNDc2MzMgMy4zNzY0NyA5LjM1NzQxIDMuNDEwNzUgOS4yNTY3OUwzLjg2MzQ3IDguMTQwOTNMMy42MTc0OSA3Ljc3NDg4TDMuNDIzNDcgNy4zNzg4M0wyLjIzMDc1IDcuMjEyOTdDMi4xMjY0NyA3LjE5MjM1IDIuMDQwNDkgNy4xMDM0MiAyLjA0MjQ1IDYuOTg2ODJMMi4wNDE4NyA1LjE4NTgyQzIuMDQzODMgNS4wNjkyMiAyLjExOTA5IDQuOTc5NTggMi4yMTcwNCA0Ljk2OTIyTDMuNDIwNjUgNC43OTM5M0wzLjg2NzQ5IDQuMDI3ODhMMy40MTEwNSAyLjkxNzMxQzMuMzczMzcgMi44MTIwNCAzLjQxMTMxIDIuNjk3NzYgMy41MTUyMyAyLjYzNzc2TDUuMDc0MDggMS43Mzc3NkM1LjE2OTM0IDEuNjgyNzYgNS4yODcyOSAxLjcwNzA0IDUuMzU5NjEgMS43OTIzMUw2LjExOTE1IDIuNzI3ODhMNi45ODAwMSAyLjczODkzTDcuNzI0OTYgMS43ODkyMkM3Ljc5MTU2IDEuNzA0NTggNy45MTU0OCAxLjY3OTIyIDguMDA4NzkgMS43NDA4Mkw5LjU2ODIxIDIuNjQxODJDOS42NzAxNyAyLjY5ODQyIDkuNzEyODUgMi44MTIzNCA5LjY4NzIzIDIuOTA3OTdMOS4yMTcxOCA0LjAzMzgzTDkuNDYzMTYgNC4zOTk4OEw5LjY1NzE4IDQuNzk1OTNaIi8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-caret-down-empty-thin: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIwIDIwIj4KCTxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSIgc2hhcGUtcmVuZGVyaW5nPSJnZW9tZXRyaWNQcmVjaXNpb24iPgoJCTxwb2x5Z29uIGNsYXNzPSJzdDEiIHBvaW50cz0iOS45LDEzLjYgMy42LDcuNCA0LjQsNi42IDkuOSwxMi4yIDE1LjQsNi43IDE2LjEsNy40ICIvPgoJPC9nPgo8L3N2Zz4K);
  --jp-icon-caret-down-empty: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDE4IDE4Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiIHNoYXBlLXJlbmRlcmluZz0iZ2VvbWV0cmljUHJlY2lzaW9uIj4KICAgIDxwYXRoIGQ9Ik01LjIsNS45TDksOS43bDMuOC0zLjhsMS4yLDEuMmwtNC45LDVsLTQuOS01TDUuMiw1Ljl6Ii8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-caret-down: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDE4IDE4Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiIHNoYXBlLXJlbmRlcmluZz0iZ2VvbWV0cmljUHJlY2lzaW9uIj4KICAgIDxwYXRoIGQ9Ik01LjIsNy41TDksMTEuMmwzLjgtMy44SDUuMnoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-caret-left: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDE4IDE4Ij4KCTxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSIgc2hhcGUtcmVuZGVyaW5nPSJnZW9tZXRyaWNQcmVjaXNpb24iPgoJCTxwYXRoIGQ9Ik0xMC44LDEyLjhMNy4xLDlsMy44LTMuOGwwLDcuNkgxMC44eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-caret-right: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDE4IDE4Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiIHNoYXBlLXJlbmRlcmluZz0iZ2VvbWV0cmljUHJlY2lzaW9uIj4KICAgIDxwYXRoIGQ9Ik03LjIsNS4yTDEwLjksOWwtMy44LDMuOFY1LjJINy4yeiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-caret-up-empty-thin: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIwIDIwIj4KCTxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSIgc2hhcGUtcmVuZGVyaW5nPSJnZW9tZXRyaWNQcmVjaXNpb24iPgoJCTxwb2x5Z29uIGNsYXNzPSJzdDEiIHBvaW50cz0iMTUuNCwxMy4zIDkuOSw3LjcgNC40LDEzLjIgMy42LDEyLjUgOS45LDYuMyAxNi4xLDEyLjYgIi8+Cgk8L2c+Cjwvc3ZnPgo=);
  --jp-icon-caret-up: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDE4IDE4Ij4KCTxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSIgc2hhcGUtcmVuZGVyaW5nPSJnZW9tZXRyaWNQcmVjaXNpb24iPgoJCTxwYXRoIGQ9Ik01LjIsMTAuNUw5LDYuOGwzLjgsMy44SDUuMnoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-case-sensitive: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIwIDIwIj4KICA8ZyBjbGFzcz0ianAtaWNvbjIiIGZpbGw9IiM0MTQxNDEiPgogICAgPHJlY3QgeD0iMiIgeT0iMiIgd2lkdGg9IjE2IiBoZWlnaHQ9IjE2Ii8+CiAgPC9nPgogIDxnIGNsYXNzPSJqcC1pY29uLWFjY2VudDIiIGZpbGw9IiNGRkYiPgogICAgPHBhdGggZD0iTTcuNiw4aDAuOWwzLjUsOGgtMS4xTDEwLDE0SDZsLTAuOSwySDRMNy42LDh6IE04LDkuMUw2LjQsMTNoMy4yTDgsOS4xeiIvPgogICAgPHBhdGggZD0iTTE2LjYsOS44Yy0wLjIsMC4xLTAuNCwwLjEtMC43LDAuMWMtMC4yLDAtMC40LTAuMS0wLjYtMC4yYy0wLjEtMC4xLTAuMi0wLjQtMC4yLTAuNyBjLTAuMywwLjMtMC42LDAuNS0wLjksMC43Yy0wLjMsMC4xLTAuNywwLjItMS4xLDAuMmMtMC4zLDAtMC41LDAtMC43LTAuMWMtMC4yLTAuMS0wLjQtMC4yLTAuNi0wLjNjLTAuMi0wLjEtMC4zLTAuMy0wLjQtMC41IGMtMC4xLTAuMi0wLjEtMC40LTAuMS0wLjdjMC0wLjMsMC4xLTAuNiwwLjItMC44YzAuMS0wLjIsMC4zLTAuNCwwLjQtMC41QzEyLDcsMTIuMiw2LjksMTIuNSw2LjhjMC4yLTAuMSwwLjUtMC4xLDAuNy0wLjIgYzAuMy0wLjEsMC41LTAuMSwwLjctMC4xYzAuMiwwLDAuNC0wLjEsMC42LTAuMWMwLjIsMCwwLjMtMC4xLDAuNC0wLjJjMC4xLTAuMSwwLjItMC4yLDAuMi0wLjRjMC0xLTEuMS0xLTEuMy0xIGMtMC40LDAtMS40LDAtMS40LDEuMmgtMC45YzAtMC40LDAuMS0wLjcsMC4yLTFjMC4xLTAuMiwwLjMtMC40LDAuNS0wLjZjMC4yLTAuMiwwLjUtMC4zLDAuOC0wLjNDMTMuMyw0LDEzLjYsNCwxMy45LDQgYzAuMywwLDAuNSwwLDAuOCwwLjFjMC4zLDAsMC41LDAuMSwwLjcsMC4yYzAuMiwwLjEsMC40LDAuMywwLjUsMC41QzE2LDUsMTYsNS4yLDE2LDUuNnYyLjljMCwwLjIsMCwwLjQsMCwwLjUgYzAsMC4xLDAuMSwwLjIsMC4zLDAuMmMwLjEsMCwwLjIsMCwwLjMsMFY5Ljh6IE0xNS4yLDYuOWMtMS4yLDAuNi0zLjEsMC4yLTMuMSwxLjRjMCwxLjQsMy4xLDEsMy4xLTAuNVY2Ljl6Ii8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-check: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTkgMTYuMTdMNC44MyAxMmwtMS40MiAxLjQxTDkgMTkgMjEgN2wtMS40MS0xLjQxeiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-circle-empty: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTEyIDJDNi40NyAyIDIgNi40NyAyIDEyczQuNDcgMTAgMTAgMTAgMTAtNC40NyAxMC0xMFMxNy41MyAyIDEyIDJ6bTAgMThjLTQuNDEgMC04LTMuNTktOC04czMuNTktOCA4LTggOCAzLjU5IDggOC0zLjU5IDgtOCA4eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-circle: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMTggMTgiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPGNpcmNsZSBjeD0iOSIgY3k9IjkiIHI9IjgiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-clear: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8bWFzayBpZD0iZG9udXRIb2xlIj4KICAgIDxyZWN0IHdpZHRoPSIyNCIgaGVpZ2h0PSIyNCIgZmlsbD0id2hpdGUiIC8+CiAgICA8Y2lyY2xlIGN4PSIxMiIgY3k9IjEyIiByPSI4IiBmaWxsPSJibGFjayIvPgogIDwvbWFzaz4KCiAgPGcgY2xhc3M9ImpwLWljb24zIiBmaWxsPSIjNjE2MTYxIj4KICAgIDxyZWN0IGhlaWdodD0iMTgiIHdpZHRoPSIyIiB4PSIxMSIgeT0iMyIgdHJhbnNmb3JtPSJyb3RhdGUoMzE1LCAxMiwgMTIpIi8+CiAgICA8Y2lyY2xlIGN4PSIxMiIgY3k9IjEyIiByPSIxMCIgbWFzaz0idXJsKCNkb251dEhvbGUpIi8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-close: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbi1ub25lIGpwLWljb24tc2VsZWN0YWJsZS1pbnZlcnNlIGpwLWljb24zLWhvdmVyIiBmaWxsPSJub25lIj4KICAgIDxjaXJjbGUgY3g9IjEyIiBjeT0iMTIiIHI9IjExIi8+CiAgPC9nPgoKICA8ZyBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIGpwLWljb24tYWNjZW50Mi1ob3ZlciIgZmlsbD0iIzYxNjE2MSI+CiAgICA8cGF0aCBkPSJNMTkgNi40MUwxNy41OSA1IDEyIDEwLjU5IDYuNDEgNSA1IDYuNDEgMTAuNTkgMTIgNSAxNy41OSA2LjQxIDE5IDEyIDEzLjQxIDE3LjU5IDE5IDE5IDE3LjU5IDEzLjQxIDEyeiIvPgogIDwvZz4KCiAgPGcgY2xhc3M9ImpwLWljb24tbm9uZSBqcC1pY29uLWJ1c3kiIGZpbGw9Im5vbmUiPgogICAgPGNpcmNsZSBjeD0iMTIiIGN5PSIxMiIgcj0iNyIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-console: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIwMCAyMDAiPgogIDxnIGNsYXNzPSJqcC1pY29uLWJyYW5kMSBqcC1pY29uLXNlbGVjdGFibGUiIGZpbGw9IiMwMjg4RDEiPgogICAgPHBhdGggZD0iTTIwIDE5LjhoMTYwdjE1OS45SDIweiIvPgogIDwvZz4KICA8ZyBjbGFzcz0ianAtaWNvbi1zZWxlY3RhYmxlLWludmVyc2UiIGZpbGw9IiNmZmYiPgogICAgPHBhdGggZD0iTTEwNSAxMjcuM2g0MHYxMi44aC00MHpNNTEuMSA3N0w3NCA5OS45bC0yMy4zIDIzLjMgMTAuNSAxMC41IDIzLjMtMjMuM0w5NSA5OS45IDg0LjUgODkuNCA2MS42IDY2LjV6Ii8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-copy: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMTggMTgiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTExLjksMUgzLjJDMi40LDEsMS43LDEuNywxLjcsMi41djEwLjJoMS41VjIuNWg4LjdWMXogTTE0LjEsMy45aC04Yy0wLjgsMC0xLjUsMC43LTEuNSwxLjV2MTAuMmMwLDAuOCwwLjcsMS41LDEuNSwxLjVoOCBjMC44LDAsMS41LTAuNywxLjUtMS41VjUuNEMxNS41LDQuNiwxNC45LDMuOSwxNC4xLDMuOXogTTE0LjEsMTUuNWgtOFY1LjRoOFYxNS41eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-cut: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTkuNjQgNy42NGMuMjMtLjUuMzYtMS4wNS4zNi0xLjY0IDAtMi4yMS0xLjc5LTQtNC00UzIgMy43OSAyIDZzMS43OSA0IDQgNGMuNTkgMCAxLjE0LS4xMyAxLjY0LS4zNkwxMCAxMmwtMi4zNiAyLjM2QzcuMTQgMTQuMTMgNi41OSAxNCA2IDE0Yy0yLjIxIDAtNCAxLjc5LTQgNHMxLjc5IDQgNCA0IDQtMS43OSA0LTRjMC0uNTktLjEzLTEuMTQtLjM2LTEuNjRMMTIgMTRsNyA3aDN2LTFMOS42NCA3LjY0ek02IDhjLTEuMSAwLTItLjg5LTItMnMuOS0yIDItMiAyIC44OSAyIDItLjkgMi0yIDJ6bTAgMTJjLTEuMSAwLTItLjg5LTItMnMuOS0yIDItMiAyIC44OSAyIDItLjkgMi0yIDJ6bTYtNy41Yy0uMjggMC0uNS0uMjItLjUtLjVzLjIyLS41LjUtLjUuNS4yMi41LjUtLjIyLjUtLjUuNXpNMTkgM2wtNiA2IDIgMiA3LTdWM3oiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-download: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTE5IDloLTRWM0g5djZINWw3IDcgNy03ek01IDE4djJoMTR2LTJINXoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-edit: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTMgMTcuMjVWMjFoMy43NUwxNy44MSA5Ljk0bC0zLjc1LTMuNzVMMyAxNy4yNXpNMjAuNzEgNy4wNGMuMzktLjM5LjM5LTEuMDIgMC0xLjQxbC0yLjM0LTIuMzRjLS4zOS0uMzktMS4wMi0uMzktMS40MSAwbC0xLjgzIDEuODMgMy43NSAzLjc1IDEuODMtMS44M3oiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-ellipses: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPGNpcmNsZSBjeD0iNSIgY3k9IjEyIiByPSIyIi8+CiAgICA8Y2lyY2xlIGN4PSIxMiIgY3k9IjEyIiByPSIyIi8+CiAgICA8Y2lyY2xlIGN4PSIxOSIgY3k9IjEyIiByPSIyIi8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-extension: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTIwLjUgMTFIMTlWN2MwLTEuMS0uOS0yLTItMmgtNFYzLjVDMTMgMi4xMiAxMS44OCAxIDEwLjUgMVM4IDIuMTIgOCAzLjVWNUg0Yy0xLjEgMC0xLjk5LjktMS45OSAydjMuOEgzLjVjMS40OSAwIDIuNyAxLjIxIDIuNyAyLjdzLTEuMjEgMi43LTIuNyAyLjdIMlYyMGMwIDEuMS45IDIgMiAyaDMuOHYtMS41YzAtMS40OSAxLjIxLTIuNyAyLjctMi43IDEuNDkgMCAyLjcgMS4yMSAyLjcgMi43VjIySDE3YzEuMSAwIDItLjkgMi0ydi00aDEuNWMxLjM4IDAgMi41LTEuMTIgMi41LTIuNVMyMS44OCAxMSAyMC41IDExeiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-fast-forward: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyNCIgaGVpZ2h0PSIyNCIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICAgIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICAgICAgPHBhdGggZD0iTTQgMThsOC41LTZMNCA2djEyem05LTEydjEybDguNS02TDEzIDZ6Ii8+CiAgICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-file-upload: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTkgMTZoNnYtNmg0bC03LTctNyA3aDR6bS00IDJoMTR2Mkg1eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-file: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMTkuMyA4LjJsLTUuNS01LjVjLS4zLS4zLS43LS41LTEuMi0uNUgzLjljLS44LjEtMS42LjktMS42IDEuOHYxNC4xYzAgLjkuNyAxLjYgMS42IDEuNmgxNC4yYy45IDAgMS42LS43IDEuNi0xLjZWOS40Yy4xLS41LS4xLS45LS40LTEuMnptLTUuOC0zLjNsMy40IDMuNmgtMy40VjQuOXptMy45IDEyLjdINC43Yy0uMSAwLS4yIDAtLjItLjJWNC43YzAtLjIuMS0uMy4yLS4zaDcuMnY0LjRzMCAuOC4zIDEuMWMuMy4zIDEuMS4zIDEuMS4zaDQuM3Y3LjJzLS4xLjItLjIuMnoiLz4KPC9zdmc+Cg==);
  --jp-icon-filter-list: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTEwIDE4aDR2LTJoLTR2MnpNMyA2djJoMThWNkgzem0zIDdoMTJ2LTJINnYyeiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-folder: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMTAgNEg0Yy0xLjEgMC0xLjk5LjktMS45OSAyTDIgMThjMCAxLjEuOSAyIDIgMmgxNmMxLjEgMCAyLS45IDItMlY4YzAtMS4xLS45LTItMi0yaC04bC0yLTJ6Ii8+Cjwvc3ZnPgo=);
  --jp-icon-html5: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDUxMiA1MTIiPgogIDxwYXRoIGNsYXNzPSJqcC1pY29uMCBqcC1pY29uLXNlbGVjdGFibGUiIGZpbGw9IiMwMDAiIGQ9Ik0xMDguNCAwaDIzdjIyLjhoMjEuMlYwaDIzdjY5aC0yM1Y0NmgtMjF2MjNoLTIzLjJNMjA2IDIzaC0yMC4zVjBoNjMuN3YyM0gyMjl2NDZoLTIzbTUzLjUtNjloMjQuMWwxNC44IDI0LjNMMzEzLjIgMGgyNC4xdjY5aC0yM1YzNC44bC0xNi4xIDI0LjgtMTYuMS0yNC44VjY5aC0yMi42bTg5LjItNjloMjN2NDYuMmgzMi42VjY5aC01NS42Ii8+CiAgPHBhdGggY2xhc3M9ImpwLWljb24tc2VsZWN0YWJsZSIgZmlsbD0iI2U0NGQyNiIgZD0iTTEwNy42IDQ3MWwtMzMtMzcwLjRoMzYyLjhsLTMzIDM3MC4yTDI1NS43IDUxMiIvPgogIDxwYXRoIGNsYXNzPSJqcC1pY29uLXNlbGVjdGFibGUiIGZpbGw9IiNmMTY1MjkiIGQ9Ik0yNTYgNDgwLjVWMTMxaDE0OC4zTDM3NiA0NDciLz4KICA8cGF0aCBjbGFzcz0ianAtaWNvbi1zZWxlY3RhYmxlLWludmVyc2UiIGZpbGw9IiNlYmViZWIiIGQ9Ik0xNDIgMTc2LjNoMTE0djQ1LjRoLTY0LjJsNC4yIDQ2LjVoNjB2NDUuM0gxNTQuNG0yIDIyLjhIMjAybDMuMiAzNi4zIDUwLjggMTMuNnY0Ny40bC05My4yLTI2Ii8+CiAgPHBhdGggY2xhc3M9ImpwLWljb24tc2VsZWN0YWJsZS1pbnZlcnNlIiBmaWxsPSIjZmZmIiBkPSJNMzY5LjYgMTc2LjNIMjU1Ljh2NDUuNGgxMDkuNm0tNC4xIDQ2LjVIMjU1Ljh2NDUuNGg1NmwtNS4zIDU5LTUwLjcgMTMuNnY0Ny4ybDkzLTI1LjgiLz4KPC9zdmc+Cg==);
  --jp-icon-image: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8cGF0aCBjbGFzcz0ianAtaWNvbi1icmFuZDQganAtaWNvbi1zZWxlY3RhYmxlLWludmVyc2UiIGZpbGw9IiNGRkYiIGQ9Ik0yLjIgMi4yaDE3LjV2MTcuNUgyLjJ6Ii8+CiAgPHBhdGggY2xhc3M9ImpwLWljb24tYnJhbmQwIGpwLWljb24tc2VsZWN0YWJsZSIgZmlsbD0iIzNGNTFCNSIgZD0iTTIuMiAyLjJ2MTcuNWgxNy41bC4xLTE3LjVIMi4yem0xMi4xIDIuMmMxLjIgMCAyLjIgMSAyLjIgMi4ycy0xIDIuMi0yLjIgMi4yLTIuMi0xLTIuMi0yLjIgMS0yLjIgMi4yLTIuMnpNNC40IDE3LjZsMy4zLTguOCAzLjMgNi42IDIuMi0zLjIgNC40IDUuNEg0LjR6Ii8+Cjwvc3ZnPgo=);
  --jp-icon-inspector: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMjAgNEg0Yy0xLjEgMC0xLjk5LjktMS45OSAyTDIgMThjMCAxLjEuOSAyIDIgMmgxNmMxLjEgMCAyLS45IDItMlY2YzAtMS4xLS45LTItMi0yem0tNSAxNEg0di00aDExdjR6bTAtNUg0VjloMTF2NHptNSA1aC00VjloNHY5eiIvPgo8L3N2Zz4K);
  --jp-icon-json: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8ZyBjbGFzcz0ianAtaWNvbi13YXJuMSBqcC1pY29uLXNlbGVjdGFibGUiIGZpbGw9IiNGOUE4MjUiPgogICAgPHBhdGggZD0iTTIwLjIgMTEuOGMtMS42IDAtMS43LjUtMS43IDEgMCAuNC4xLjkuMSAxLjMuMS41LjEuOS4xIDEuMyAwIDEuNy0xLjQgMi4zLTMuNSAyLjNoLS45di0xLjloLjVjMS4xIDAgMS40IDAgMS40LS44IDAtLjMgMC0uNi0uMS0xIDAtLjQtLjEtLjgtLjEtMS4yIDAtMS4zIDAtMS44IDEuMy0yLTEuMy0uMi0xLjMtLjctMS4zLTIgMC0uNC4xLS44LjEtMS4yLjEtLjQuMS0uNy4xLTEgMC0uOC0uNC0uNy0xLjQtLjhoLS41VjQuMWguOWMyLjIgMCAzLjUuNyAzLjUgMi4zIDAgLjQtLjEuOS0uMSAxLjMtLjEuNS0uMS45LS4xIDEuMyAwIC41LjIgMSAxLjcgMXYxLjh6TTEuOCAxMC4xYzEuNiAwIDEuNy0uNSAxLjctMSAwLS40LS4xLS45LS4xLTEuMy0uMS0uNS0uMS0uOS0uMS0xLjMgMC0xLjYgMS40LTIuMyAzLjUtMi4zaC45djEuOWgtLjVjLTEgMC0xLjQgMC0xLjQuOCAwIC4zIDAgLjYuMSAxIDAgLjIuMS42LjEgMSAwIDEuMyAwIDEuOC0xLjMgMkM2IDExLjIgNiAxMS43IDYgMTNjMCAuNC0uMS44LS4xIDEuMi0uMS4zLS4xLjctLjEgMSAwIC44LjMuOCAxLjQuOGguNXYxLjloLS45Yy0yLjEgMC0zLjUtLjYtMy41LTIuMyAwLS40LjEtLjkuMS0xLjMuMS0uNS4xLS45LjEtMS4zIDAtLjUtLjItMS0xLjctMXYtMS45eiIvPgogICAgPGNpcmNsZSBjeD0iMTEiIGN5PSIxMy44IiByPSIyLjEiLz4KICAgIDxjaXJjbGUgY3g9IjExIiBjeT0iOC4yIiByPSIyLjEiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-jupyter-favicon: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTUyIiBoZWlnaHQ9IjE2NSIgdmlld0JveD0iMCAwIDE1MiAxNjUiIHZlcnNpb249IjEuMSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbi13YXJuMCIgZmlsbD0iI0YzNzcyNiI+CiAgICA8cGF0aCB0cmFuc2Zvcm09InRyYW5zbGF0ZSgwLjA3ODk0NywgMTEwLjU4MjkyNykiIGQ9Ik03NS45NDIyODQyLDI5LjU4MDQ1NjEgQzQzLjMwMjM5NDcsMjkuNTgwNDU2MSAxNC43OTY3ODMyLDE3LjY1MzQ2MzQgMCwwIEM1LjUxMDgzMjExLDE1Ljg0MDY4MjkgMTUuNzgxNTM4OSwyOS41NjY3NzMyIDI5LjM5MDQ5NDcsMzkuMjc4NDE3MSBDNDIuOTk5Nyw0OC45ODk4NTM3IDU5LjI3MzcsNTQuMjA2NzgwNSA3NS45NjA1Nzg5LDU0LjIwNjc4MDUgQzkyLjY0NzQ1NzksNTQuMjA2NzgwNSAxMDguOTIxNDU4LDQ4Ljk4OTg1MzcgMTIyLjUzMDY2MywzOS4yNzg0MTcxIEMxMzYuMTM5NDUzLDI5LjU2Njc3MzIgMTQ2LjQxMDI4NCwxNS44NDA2ODI5IDE1MS45MjExNTgsMCBDMTM3LjA4Nzg2OCwxNy42NTM0NjM0IDEwOC41ODI1ODksMjkuNTgwNDU2MSA3NS45NDIyODQyLDI5LjU4MDQ1NjEgTDc1Ljk0MjI4NDIsMjkuNTgwNDU2MSBaIiAvPgogICAgPHBhdGggdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMC4wMzczNjgsIDAuNzA0ODc4KSIgZD0iTTc1Ljk3ODQ1NzksMjQuNjI2NDA3MyBDMTA4LjYxODc2MywyNC42MjY0MDczIDEzNy4xMjQ0NTgsMzYuNTUzNDQxNSAxNTEuOTIxMTU4LDU0LjIwNjc4MDUgQzE0Ni40MTAyODQsMzguMzY2MjIyIDEzNi4xMzk0NTMsMjQuNjQwMTMxNyAxMjIuNTMwNjYzLDE0LjkyODQ4NzggQzEwOC45MjE0NTgsNS4yMTY4NDM5IDkyLjY0NzQ1NzksMCA3NS45NjA1Nzg5LDAgQzU5LjI3MzcsMCA0Mi45OTk3LDUuMjE2ODQzOSAyOS4zOTA0OTQ3LDE0LjkyODQ4NzggQzE1Ljc4MTUzODksMjQuNjQwMTMxNyA1LjUxMDgzMjExLDM4LjM2NjIyMiAwLDU0LjIwNjc4MDUgQzE0LjgzMzA4MTYsMzYuNTg5OTI5MyA0My4zMzg1Njg0LDI0LjYyNjQwNzMgNzUuOTc4NDU3OSwyNC42MjY0MDczIEw3NS45Nzg0NTc5LDI0LjYyNjQwNzMgWiIgLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-jupyter: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzkiIGhlaWdodD0iNTEiIHZpZXdCb3g9IjAgMCAzOSA1MSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMTYzOCAtMjI4MSkiPgogICAgPGcgY2xhc3M9ImpwLWljb24td2FybjAiIGZpbGw9IiNGMzc3MjYiPgogICAgICA8cGF0aCB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxNjM5Ljc0IDIzMTEuOTgpIiBkPSJNIDE4LjI2NDYgNy4xMzQxMUMgMTAuNDE0NSA3LjEzNDExIDMuNTU4NzIgNC4yNTc2IDAgMEMgMS4zMjUzOSAzLjgyMDQgMy43OTU1NiA3LjEzMDgxIDcuMDY4NiA5LjQ3MzAzQyAxMC4zNDE3IDExLjgxNTIgMTQuMjU1NyAxMy4wNzM0IDE4LjI2OSAxMy4wNzM0QyAyMi4yODIzIDEzLjA3MzQgMjYuMTk2MyAxMS44MTUyIDI5LjQ2OTQgOS40NzMwM0MgMzIuNzQyNCA3LjEzMDgxIDM1LjIxMjYgMy44MjA0IDM2LjUzOCAwQyAzMi45NzA1IDQuMjU3NiAyNi4xMTQ4IDcuMTM0MTEgMTguMjY0NiA3LjEzNDExWiIvPgogICAgICA8cGF0aCB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxNjM5LjczIDIyODUuNDgpIiBkPSJNIDE4LjI3MzMgNS45MzkzMUMgMjYuMTIzNSA1LjkzOTMxIDMyLjk3OTMgOC44MTU4MyAzNi41MzggMTMuMDczNEMgMzUuMjEyNiA5LjI1MzAzIDMyLjc0MjQgNS45NDI2MiAyOS40Njk0IDMuNjAwNEMgMjYuMTk2MyAxLjI1ODE4IDIyLjI4MjMgMCAxOC4yNjkgMEMgMTQuMjU1NyAwIDEwLjM0MTcgMS4yNTgxOCA3LjA2ODYgMy42MDA0QyAzLjc5NTU2IDUuOTQyNjIgMS4zMjUzOSA5LjI1MzAzIDAgMTMuMDczNEMgMy41Njc0NSA4LjgyNDYzIDEwLjQyMzIgNS45MzkzMSAxOC4yNzMzIDUuOTM5MzFaIi8+CiAgICA8L2c+CiAgICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgICA8cGF0aCB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxNjY5LjMgMjI4MS4zMSkiIGQ9Ik0gNS44OTM1MyAyLjg0NEMgNS45MTg4OSAzLjQzMTY1IDUuNzcwODUgNC4wMTM2NyA1LjQ2ODE1IDQuNTE2NDVDIDUuMTY1NDUgNS4wMTkyMiA0LjcyMTY4IDUuNDIwMTUgNC4xOTI5OSA1LjY2ODUxQyAzLjY2NDMgNS45MTY4OCAzLjA3NDQ0IDYuMDAxNTEgMi40OTgwNSA1LjkxMTcxQyAxLjkyMTY2IDUuODIxOSAxLjM4NDYzIDUuNTYxNyAwLjk1NDg5OCA1LjE2NDAxQyAwLjUyNTE3IDQuNzY2MzMgMC4yMjIwNTYgNC4yNDkwMyAwLjA4MzkwMzcgMy42Nzc1N0MgLTAuMDU0MjQ4MyAzLjEwNjExIC0wLjAyMTIzIDIuNTA2MTcgMC4xNzg3ODEgMS45NTM2NEMgMC4zNzg3OTMgMS40MDExIDAuNzM2ODA5IDAuOTIwODE3IDEuMjA3NTQgMC41NzM1MzhDIDEuNjc4MjYgMC4yMjYyNTkgMi4yNDA1NSAwLjAyNzU5MTkgMi44MjMyNiAwLjAwMjY3MjI5QyAzLjYwMzg5IC0wLjAzMDcxMTUgNC4zNjU3MyAwLjI0OTc4OSA0Ljk0MTQyIDAuNzgyNTUxQyA1LjUxNzExIDEuMzE1MzEgNS44NTk1NiAyLjA1Njc2IDUuODkzNTMgMi44NDRaIi8+CiAgICAgIDxwYXRoIHRyYW5zZm9ybT0idHJhbnNsYXRlKDE2MzkuOCAyMzIzLjgxKSIgZD0iTSA3LjQyNzg5IDMuNTgzMzhDIDcuNDYwMDggNC4zMjQzIDcuMjczNTUgNS4wNTgxOSA2Ljg5MTkzIDUuNjkyMTNDIDYuNTEwMzEgNi4zMjYwNyA1Ljk1MDc1IDYuODMxNTYgNS4yODQxMSA3LjE0NDZDIDQuNjE3NDcgNy40NTc2MyAzLjg3MzcxIDcuNTY0MTQgMy4xNDcwMiA3LjQ1MDYzQyAyLjQyMDMyIDcuMzM3MTIgMS43NDMzNiA3LjAwODcgMS4yMDE4NCA2LjUwNjk1QyAwLjY2MDMyOCA2LjAwNTIgMC4yNzg2MSA1LjM1MjY4IDAuMTA1MDE3IDQuNjMyMDJDIC0wLjA2ODU3NTcgMy45MTEzNSAtMC4wMjYyMzYxIDMuMTU0OTQgMC4yMjY2NzUgMi40NTg1NkMgMC40Nzk1ODcgMS43NjIxNyAwLjkzMTY5NyAxLjE1NzEzIDEuNTI1NzYgMC43MjAwMzNDIDIuMTE5ODMgMC4yODI5MzUgMi44MjkxNCAwLjAzMzQzOTUgMy41NjM4OSAwLjAwMzEzMzQ0QyA0LjU0NjY3IC0wLjAzNzQwMzMgNS41MDUyOSAwLjMxNjcwNiA2LjIyOTYxIDAuOTg3ODM1QyA2Ljk1MzkzIDEuNjU4OTYgNy4zODQ4NCAyLjU5MjM1IDcuNDI3ODkgMy41ODMzOEwgNy40Mjc4OSAzLjU4MzM4WiIvPgogICAgICA8cGF0aCB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxNjM4LjM2IDIyODYuMDYpIiBkPSJNIDIuMjc0NzEgNC4zOTYyOUMgMS44NDM2MyA0LjQxNTA4IDEuNDE2NzEgNC4zMDQ0NSAxLjA0Nzk5IDQuMDc4NDNDIDAuNjc5MjY4IDMuODUyNCAwLjM4NTMyOCAzLjUyMTE0IDAuMjAzMzcxIDMuMTI2NTZDIDAuMDIxNDEzNiAyLjczMTk4IC0wLjA0MDM3OTggMi4yOTE4MyAwLjAyNTgxMTYgMS44NjE4MUMgMC4wOTIwMDMxIDEuNDMxOCAwLjI4MzIwNCAxLjAzMTI2IDAuNTc1MjEzIDAuNzEwODgzQyAwLjg2NzIyMiAwLjM5MDUxIDEuMjQ2OTEgMC4xNjQ3MDggMS42NjYyMiAwLjA2MjA1OTJDIDIuMDg1NTMgLTAuMDQwNTg5NyAyLjUyNTYxIC0wLjAxNTQ3MTQgMi45MzA3NiAwLjEzNDIzNUMgMy4zMzU5MSAwLjI4Mzk0MSAzLjY4NzkyIDAuNTUxNTA1IDMuOTQyMjIgMC45MDMwNkMgNC4xOTY1MiAxLjI1NDYyIDQuMzQxNjkgMS42NzQzNiA0LjM1OTM1IDIuMTA5MTZDIDQuMzgyOTkgMi42OTEwNyA0LjE3Njc4IDMuMjU4NjkgMy43ODU5NyAzLjY4NzQ2QyAzLjM5NTE2IDQuMTE2MjQgMi44NTE2NiA0LjM3MTE2IDIuMjc0NzEgNC4zOTYyOUwgMi4yNzQ3MSA0LjM5NjI5WiIvPgogICAgPC9nPgogIDwvZz4+Cjwvc3ZnPgo=);
  --jp-icon-jupyterlab-wordmark: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyMDAiIHZpZXdCb3g9IjAgMCAxODYwLjggNDc1Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjIiIGZpbGw9IiM0RTRFNEUiIHRyYW5zZm9ybT0idHJhbnNsYXRlKDQ4MC4xMzY0MDEsIDY0LjI3MTQ5MykiPgogICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMC4wMDAwMDAsIDU4Ljg3NTU2NikiPgogICAgICA8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgwLjA4NzYwMywgMC4xNDAyOTQpIj4KICAgICAgICA8cGF0aCBkPSJNLTQyNi45LDE2OS44YzAsNDguNy0zLjcsNjQuNy0xMy42LDc2LjRjLTEwLjgsMTAtMjUsMTUuNS0zOS43LDE1LjVsMy43LDI5IGMyMi44LDAuMyw0NC44LTcuOSw2MS45LTIzLjFjMTcuOC0xOC41LDI0LTQ0LjEsMjQtODMuM1YwSC00Mjd2MTcwLjFMLTQyNi45LDE2OS44TC00MjYuOSwxNjkuOHoiLz4KICAgICAgPC9nPgogICAgPC9nPgogICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMTU1LjA0NTI5NiwgNTYuODM3MTA0KSI+CiAgICAgIDxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDEuNTYyNDUzLCAxLjc5OTg0MikiPgogICAgICAgIDxwYXRoIGQ9Ik0tMzEyLDE0OGMwLDIxLDAsMzkuNSwxLjcsNTUuNGgtMzEuOGwtMi4xLTMzLjNoLTAuOGMtNi43LDExLjYtMTYuNCwyMS4zLTI4LDI3LjkgYy0xMS42LDYuNi0yNC44LDEwLTM4LjIsOS44Yy0zMS40LDAtNjktMTcuNy02OS04OVYwaDM2LjR2MTEyLjdjMCwzOC43LDExLjYsNjQuNyw0NC42LDY0LjdjMTAuMy0wLjIsMjAuNC0zLjUsMjguOS05LjQgYzguNS01LjksMTUuMS0xNC4zLDE4LjktMjMuOWMyLjItNi4xLDMuMy0xMi41LDMuMy0xOC45VjAuMmgzNi40VjE0OEgtMzEyTC0zMTIsMTQ4eiIvPgogICAgICA8L2c+CiAgICA8L2c+CiAgICA8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgzOTAuMDEzMzIyLCA1My40Nzk2MzgpIj4KICAgICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMS43MDY0NTgsIDAuMjMxNDI1KSI+CiAgICAgICAgPHBhdGggZD0iTS00NzguNiw3MS40YzAtMjYtMC44LTQ3LTEuNy02Ni43aDMyLjdsMS43LDM0LjhoMC44YzcuMS0xMi41LDE3LjUtMjIuOCwzMC4xLTI5LjcgYzEyLjUtNywyNi43LTEwLjMsNDEtOS44YzQ4LjMsMCw4NC43LDQxLjcsODQuNywxMDMuM2MwLDczLjEtNDMuNywxMDkuMi05MSwxMDkuMmMtMTIuMSwwLjUtMjQuMi0yLjItMzUtNy44IGMtMTAuOC01LjYtMTkuOS0xMy45LTI2LjYtMjQuMmgtMC44VjI5MWgtMzZ2LTIyMEwtNDc4LjYsNzEuNEwtNDc4LjYsNzEuNHogTS00NDIuNiwxMjUuNmMwLjEsNS4xLDAuNiwxMC4xLDEuNywxNS4xIGMzLDEyLjMsOS45LDIzLjMsMTkuOCwzMS4xYzkuOSw3LjgsMjIuMSwxMi4xLDM0LjcsMTIuMWMzOC41LDAsNjAuNy0zMS45LDYwLjctNzguNWMwLTQwLjctMjEuMS03NS42LTU5LjUtNzUuNiBjLTEyLjksMC40LTI1LjMsNS4xLTM1LjMsMTMuNGMtOS45LDguMy0xNi45LDE5LjctMTkuNiwzMi40Yy0xLjUsNC45LTIuMywxMC0yLjUsMTUuMVYxMjUuNkwtNDQyLjYsMTI1LjZMLTQ0Mi42LDEyNS42eiIvPgogICAgICA8L2c+CiAgICA8L2c+CiAgICA8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSg2MDYuNzQwNzI2LCA1Ni44MzcxMDQpIj4KICAgICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMC43NTEyMjYsIDEuOTg5Mjk5KSI+CiAgICAgICAgPHBhdGggZD0iTS00NDAuOCwwbDQzLjcsMTIwLjFjNC41LDEzLjQsOS41LDI5LjQsMTIuOCw0MS43aDAuOGMzLjctMTIuMiw3LjktMjcuNywxMi44LTQyLjQgbDM5LjctMTE5LjJoMzguNUwtMzQ2LjksMTQ1Yy0yNiw2OS43LTQzLjcsMTA1LjQtNjguNiwxMjcuMmMtMTIuNSwxMS43LTI3LjksMjAtNDQuNiwyMy45bC05LjEtMzEuMSBjMTEuNy0zLjksMjIuNS0xMC4xLDMxLjgtMTguMWMxMy4yLTExLjEsMjMuNy0yNS4yLDMwLjYtNDEuMmMxLjUtMi44LDIuNS01LjcsMi45LTguOGMtMC4zLTMuMy0xLjItNi42LTIuNS05LjdMLTQ4MC4yLDAuMSBoMzkuN0wtNDQwLjgsMEwtNDQwLjgsMHoiLz4KICAgICAgPC9nPgogICAgPC9nPgogICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoODIyLjc0ODEwNCwgMC4wMDAwMDApIj4KICAgICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMS40NjQwNTAsIDAuMzc4OTE0KSI+CiAgICAgICAgPHBhdGggZD0iTS00MTMuNywwdjU4LjNoNTJ2MjguMmgtNTJWMTk2YzAsMjUsNywzOS41LDI3LjMsMzkuNWM3LjEsMC4xLDE0LjItMC43LDIxLjEtMi41IGwxLjcsMjcuN2MtMTAuMywzLjctMjEuMyw1LjQtMzIuMiw1Yy03LjMsMC40LTE0LjYtMC43LTIxLjMtMy40Yy02LjgtMi43LTEyLjktNi44LTE3LjktMTIuMWMtMTAuMy0xMC45LTE0LjEtMjktMTQuMS01Mi45IFY4Ni41aC0zMVY1OC4zaDMxVjkuNkwtNDEzLjcsMEwtNDEzLjcsMHoiLz4KICAgICAgPC9nPgogICAgPC9nPgogICAgPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoOTc0LjQzMzI4NiwgNTMuNDc5NjM4KSI+CiAgICAgIDxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAuOTkwMDM0LCAwLjYxMDMzOSkiPgogICAgICAgIDxwYXRoIGQ9Ik0tNDQ1LjgsMTEzYzAuOCw1MCwzMi4yLDcwLjYsNjguNiw3MC42YzE5LDAuNiwzNy45LTMsNTUuMy0xMC41bDYuMiwyNi40IGMtMjAuOSw4LjktNDMuNSwxMy4xLTY2LjIsMTIuNmMtNjEuNSwwLTk4LjMtNDEuMi05OC4zLTEwMi41Qy00ODAuMiw0OC4yLTQ0NC43LDAtMzg2LjUsMGM2NS4yLDAsODIuNyw1OC4zLDgyLjcsOTUuNyBjLTAuMSw1LjgtMC41LDExLjUtMS4yLDE3LjJoLTE0MC42SC00NDUuOEwtNDQ1LjgsMTEzeiBNLTMzOS4yLDg2LjZjMC40LTIzLjUtOS41LTYwLjEtNTAuNC02MC4xIGMtMzYuOCwwLTUyLjgsMzQuNC01NS43LDYwLjFILTMzOS4yTC0zMzkuMiw4Ni42TC0zMzkuMiw4Ni42eiIvPgogICAgICA8L2c+CiAgICA8L2c+CiAgICA8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxMjAxLjk2MTA1OCwgNTMuNDc5NjM4KSI+CiAgICAgIDxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDEuMTc5NjQwLCAwLjcwNTA2OCkiPgogICAgICAgIDxwYXRoIGQ9Ik0tNDc4LjYsNjhjMC0yMy45LTAuNC00NC41LTEuNy02My40aDMxLjhsMS4yLDM5LjloMS43YzkuMS0yNy4zLDMxLTQ0LjUsNTUuMy00NC41IGMzLjUtMC4xLDcsMC40LDEwLjMsMS4ydjM0LjhjLTQuMS0wLjktOC4yLTEuMy0xMi40LTEuMmMtMjUuNiwwLTQzLjcsMTkuNy00OC43LDQ3LjRjLTEsNS43LTEuNiwxMS41LTEuNywxNy4ydjEwOC4zaC0zNlY2OCBMLTQ3OC42LDY4eiIvPgogICAgICA8L2c+CiAgICA8L2c+CiAgPC9nPgoKICA8ZyBjbGFzcz0ianAtaWNvbi13YXJuMCIgZmlsbD0iI0YzNzcyNiI+CiAgICA8cGF0aCBkPSJNMTM1Mi4zLDMyNi4yaDM3VjI4aC0zN1YzMjYuMnogTTE2MDQuOCwzMjYuMmMtMi41LTEzLjktMy40LTMxLjEtMy40LTQ4Ljd2LTc2IGMwLTQwLjctMTUuMS04My4xLTc3LjMtODMuMWMtMjUuNiwwLTUwLDcuMS02Ni44LDE4LjFsOC40LDI0LjRjMTQuMy05LjIsMzQtMTUuMSw1My0xNS4xYzQxLjYsMCw0Ni4yLDMwLjIsNDYuMiw0N3Y0LjIgYy03OC42LTAuNC0xMjIuMywyNi41LTEyMi4zLDc1LjZjMCwyOS40LDIxLDU4LjQsNjIuMiw1OC40YzI5LDAsNTAuOS0xNC4zLDYyLjItMzAuMmgxLjNsMi45LDI1LjZIMTYwNC44eiBNMTU2NS43LDI1Ny43IGMwLDMuOC0wLjgsOC0yLjEsMTEuOGMtNS45LDE3LjItMjIuNywzNC00OS4yLDM0Yy0xOC45LDAtMzQuOS0xMS4zLTM0LjktMzUuM2MwLTM5LjUsNDUuOC00Ni42LDg2LjItNDUuOFYyNTcuN3ogTTE2OTguNSwzMjYuMiBsMS43LTMzLjZoMS4zYzE1LjEsMjYuOSwzOC43LDM4LjIsNjguMSwzOC4yYzQ1LjQsMCw5MS4yLTM2LjEsOTEuMi0xMDguOGMwLjQtNjEuNy0zNS4zLTEwMy43LTg1LjctMTAzLjcgYy0zMi44LDAtNTYuMywxNC43LTY5LjMsMzcuNGgtMC44VjI4aC0zNi42djI0NS43YzAsMTguMS0wLjgsMzguNi0xLjcsNTIuNUgxNjk4LjV6IE0xNzA0LjgsMjA4LjJjMC01LjksMS4zLTEwLjksMi4xLTE1LjEgYzcuNi0yOC4xLDMxLjEtNDUuNCw1Ni4zLTQ1LjRjMzkuNSwwLDYwLjUsMzQuOSw2MC41LDc1LjZjMCw0Ni42LTIzLjEsNzguMS02MS44LDc4LjFjLTI2LjksMC00OC4zLTE3LjYtNTUuNS00My4zIGMtMC44LTQuMi0xLjctOC44LTEuNy0xMy40VjIwOC4yeiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-kernel: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICAgIDxwYXRoIGNsYXNzPSJqcC1pY29uMiIgZmlsbD0iIzYxNjE2MSIgZD0iTTE1IDlIOXY2aDZWOXptLTIgNGgtMnYtMmgydjJ6bTgtMlY5aC0yVjdjMC0xLjEtLjktMi0yLTJoLTJWM2gtMnYyaC0yVjNIOXYySDdjLTEuMSAwLTIgLjktMiAydjJIM3YyaDJ2MkgzdjJoMnYyYzAgMS4xLjkgMiAyIDJoMnYyaDJ2LTJoMnYyaDJ2LTJoMmMxLjEgMCAyLS45IDItMnYtMmgydi0yaC0ydi0yaDJ6bS00IDZIN1Y3aDEwdjEweiIvPgo8L3N2Zz4K);
  --jp-icon-keyboard: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMjAgNUg0Yy0xLjEgMC0xLjk5LjktMS45OSAyTDIgMTdjMCAxLjEuOSAyIDIgMmgxNmMxLjEgMCAyLS45IDItMlY3YzAtMS4xLS45LTItMi0yem0tOSAzaDJ2MmgtMlY4em0wIDNoMnYyaC0ydi0yek04IDhoMnYySDhWOHptMCAzaDJ2Mkg4di0yem0tMSAySDV2LTJoMnYyem0wLTNINVY4aDJ2MnptOSA3SDh2LTJoOHYyem0wLTRoLTJ2LTJoMnYyem0wLTNoLTJWOGgydjJ6bTMgM2gtMnYtMmgydjJ6bTAtM2gtMlY4aDJ2MnoiLz4KPC9zdmc+Cg==);
  --jp-icon-launcher: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMTkgMTlINVY1aDdWM0g1YTIgMiAwIDAwLTIgMnYxNGEyIDIgMCAwMDIgMmgxNGMxLjEgMCAyLS45IDItMnYtN2gtMnY3ek0xNCAzdjJoMy41OWwtOS44MyA5LjgzIDEuNDEgMS40MUwxOSA2LjQxVjEwaDJWM2gtN3oiLz4KPC9zdmc+Cg==);
  --jp-icon-line-form: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICAgIDxwYXRoIGZpbGw9IndoaXRlIiBkPSJNNS44OCA0LjEyTDEzLjc2IDEybC03Ljg4IDcuODhMOCAyMmwxMC0xMEw4IDJ6Ii8+Cjwvc3ZnPgo=);
  --jp-icon-link: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTMuOSAxMmMwLTEuNzEgMS4zOS0zLjEgMy4xLTMuMWg0VjdIN2MtMi43NiAwLTUgMi4yNC01IDVzMi4yNCA1IDUgNWg0di0xLjlIN2MtMS43MSAwLTMuMS0xLjM5LTMuMS0zLjF6TTggMTNoOHYtMkg4djJ6bTktNmgtNHYxLjloNGMxLjcxIDAgMy4xIDEuMzkgMy4xIDMuMXMtMS4zOSAzLjEtMy4xIDMuMWgtNFYxN2g0YzIuNzYgMCA1LTIuMjQgNS01cy0yLjI0LTUtNS01eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-list: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICAgIDxwYXRoIGNsYXNzPSJqcC1pY29uMiBqcC1pY29uLXNlbGVjdGFibGUiIGZpbGw9IiM2MTYxNjEiIGQ9Ik0xOSA1djE0SDVWNWgxNG0xLjEtMkgzLjljLS41IDAtLjkuNC0uOS45djE2LjJjMCAuNC40LjkuOS45aDE2LjJjLjQgMCAuOS0uNS45LS45VjMuOWMwLS41LS41LS45LS45LS45ek0xMSA3aDZ2MmgtNlY3em0wIDRoNnYyaC02di0yem0wIDRoNnYyaC02ek03IDdoMnYySDd6bTAgNGgydjJIN3ptMCA0aDJ2Mkg3eiIvPgo8L3N2Zz4=);
  --jp-icon-listings-info: url(data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iaXNvLTg4NTktMSI/Pg0KPHN2ZyB2ZXJzaW9uPSIxLjEiIGlkPSJDYXBhXzEiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHg9IjBweCIgeT0iMHB4Ig0KCSB2aWV3Qm94PSIwIDAgNTAuOTc4IDUwLjk3OCIgc3R5bGU9ImVuYWJsZS1iYWNrZ3JvdW5kOm5ldyAwIDAgNTAuOTc4IDUwLjk3ODsiIHhtbDpzcGFjZT0icHJlc2VydmUiPg0KPGc+DQoJPGc+DQoJCTxnPg0KCQkJPHBhdGggc3R5bGU9ImZpbGw6IzAxMDAwMjsiIGQ9Ik00My41Miw3LjQ1OEMzOC43MTEsMi42NDgsMzIuMzA3LDAsMjUuNDg5LDBDMTguNjcsMCwxMi4yNjYsMi42NDgsNy40NTgsNy40NTgNCgkJCQljLTkuOTQzLDkuOTQxLTkuOTQzLDI2LjExOSwwLDM2LjA2MmM0LjgwOSw0LjgwOSwxMS4yMTIsNy40NTYsMTguMDMxLDcuNDU4YzAsMCwwLjAwMSwwLDAuMDAyLDANCgkJCQljNi44MTYsMCwxMy4yMjEtMi42NDgsMTguMDI5LTcuNDU4YzQuODA5LTQuODA5LDcuNDU3LTExLjIxMiw3LjQ1Ny0xOC4wM0M1MC45NzcsMTguNjcsNDguMzI4LDEyLjI2Niw0My41Miw3LjQ1OHoNCgkJCQkgTTQyLjEwNiw0Mi4xMDVjLTQuNDMyLDQuNDMxLTEwLjMzMiw2Ljg3Mi0xNi42MTUsNi44NzJoLTAuMDAyYy02LjI4NS0wLjAwMS0xMi4xODctMi40NDEtMTYuNjE3LTYuODcyDQoJCQkJYy05LjE2Mi05LjE2My05LjE2Mi0yNC4wNzEsMC0zMy4yMzNDMTMuMzAzLDQuNDQsMTkuMjA0LDIsMjUuNDg5LDJjNi4yODQsMCwxMi4xODYsMi40NCwxNi42MTcsNi44NzINCgkJCQljNC40MzEsNC40MzEsNi44NzEsMTAuMzMyLDYuODcxLDE2LjYxN0M0OC45NzcsMzEuNzcyLDQ2LjUzNiwzNy42NzUsNDIuMTA2LDQyLjEwNXoiLz4NCgkJPC9nPg0KCQk8Zz4NCgkJCTxwYXRoIHN0eWxlPSJmaWxsOiMwMTAwMDI7IiBkPSJNMjMuNTc4LDMyLjIxOGMtMC4wMjMtMS43MzQsMC4xNDMtMy4wNTksMC40OTYtMy45NzJjMC4zNTMtMC45MTMsMS4xMS0xLjk5NywyLjI3Mi0zLjI1Mw0KCQkJCWMwLjQ2OC0wLjUzNiwwLjkyMy0xLjA2MiwxLjM2Ny0xLjU3NWMwLjYyNi0wLjc1MywxLjEwNC0xLjQ3OCwxLjQzNi0yLjE3NWMwLjMzMS0wLjcwNywwLjQ5NS0xLjU0MSwwLjQ5NS0yLjUNCgkJCQljMC0xLjA5Ni0wLjI2LTIuMDg4LTAuNzc5LTIuOTc5Yy0wLjU2NS0wLjg3OS0xLjUwMS0xLjMzNi0yLjgwNi0xLjM2OWMtMS44MDIsMC4wNTctMi45ODUsMC42NjctMy41NSwxLjgzMg0KCQkJCWMtMC4zMDEsMC41MzUtMC41MDMsMS4xNDEtMC42MDcsMS44MTRjLTAuMTM5LDAuNzA3LTAuMjA3LDEuNDMyLTAuMjA3LDIuMTc0aC0yLjkzN2MtMC4wOTEtMi4yMDgsMC40MDctNC4xMTQsMS40OTMtNS43MTkNCgkJCQljMS4wNjItMS42NCwyLjg1NS0yLjQ4MSw1LjM3OC0yLjUyN2MyLjE2LDAuMDIzLDMuODc0LDAuNjA4LDUuMTQxLDEuNzU4YzEuMjc4LDEuMTYsMS45MjksMi43NjQsMS45NSw0LjgxMQ0KCQkJCWMwLDEuMTQyLTAuMTM3LDIuMTExLTAuNDEsMi45MTFjLTAuMzA5LDAuODQ1LTAuNzMxLDEuNTkzLTEuMjY4LDIuMjQzYy0wLjQ5MiwwLjY1LTEuMDY4LDEuMzE4LTEuNzMsMi4wMDINCgkJCQljLTAuNjUsMC42OTctMS4zMTMsMS40NzktMS45ODcsMi4zNDZjLTAuMjM5LDAuMzc3LTAuNDI5LDAuNzc3LTAuNTY1LDEuMTk5Yy0wLjE2LDAuOTU5LTAuMjE3LDEuOTUxLTAuMTcxLDIuOTc5DQoJCQkJQzI2LjU4OSwzMi4yMTgsMjMuNTc4LDMyLjIxOCwyMy41NzgsMzIuMjE4eiBNMjMuNTc4LDM4LjIydi0zLjQ4NGgzLjA3NnYzLjQ4NEgyMy41Nzh6Ii8+DQoJCTwvZz4NCgk8L2c+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8L3N2Zz4NCg==);
  --jp-icon-markdown: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8cGF0aCBjbGFzcz0ianAtaWNvbi1jb250cmFzdDAganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjN0IxRkEyIiBkPSJNNSAxNC45aDEybC02LjEgNnptOS40LTYuOGMwLTEuMy0uMS0yLjktLjEtNC41LS40IDEuNC0uOSAyLjktMS4zIDQuM2wtMS4zIDQuM2gtMkw4LjUgNy45Yy0uNC0xLjMtLjctMi45LTEtNC4zLS4xIDEuNi0uMSAzLjItLjIgNC42TDcgMTIuNEg0LjhsLjctMTFoMy4zTDEwIDVjLjQgMS4yLjcgMi43IDEgMy45LjMtMS4yLjctMi42IDEtMy45bDEuMi0zLjdoMy4zbC42IDExaC0yLjRsLS4zLTQuMnoiLz4KPC9zdmc+Cg==);
  --jp-icon-new-folder: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTIwIDZoLThsLTItMkg0Yy0xLjExIDAtMS45OS44OS0xLjk5IDJMMiAxOGMwIDEuMTEuODkgMiAyIDJoMTZjMS4xMSAwIDItLjg5IDItMlY4YzAtMS4xMS0uODktMi0yLTJ6bS0xIDhoLTN2M2gtMnYtM2gtM3YtMmgzVjloMnYzaDN2MnoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-not-trusted: url(data:image/svg+xml;base64,PHN2ZyBmaWxsPSJub25lIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI1IDI1Ij4KICAgIDxwYXRoIGNsYXNzPSJqcC1pY29uMiIgc3Ryb2tlPSIjMzMzMzMzIiBzdHJva2Utd2lkdGg9IjIiIHRyYW5zZm9ybT0idHJhbnNsYXRlKDMgMykiIGQ9Ik0xLjg2MDk0IDExLjQ0MDlDMC44MjY0NDggOC43NzAyNyAwLjg2Mzc3OSA2LjA1NzY0IDEuMjQ5MDcgNC4xOTkzMkMyLjQ4MjA2IDMuOTMzNDcgNC4wODA2OCAzLjQwMzQ3IDUuNjAxMDIgMi44NDQ5QzcuMjM1NDkgMi4yNDQ0IDguODU2NjYgMS41ODE1IDkuOTg3NiAxLjA5NTM5QzExLjA1OTcgMS41ODM0MSAxMi42MDk0IDIuMjQ0NCAxNC4yMTggMi44NDMzOUMxNS43NTAzIDMuNDEzOTQgMTcuMzk5NSAzLjk1MjU4IDE4Ljc1MzkgNC4yMTM4NUMxOS4xMzY0IDYuMDcxNzcgMTkuMTcwOSA4Ljc3NzIyIDE4LjEzOSAxMS40NDA5QzE3LjAzMDMgMTQuMzAzMiAxNC42NjY4IDE3LjE4NDQgOS45OTk5OSAxOC45MzU0QzUuMzMzMTkgMTcuMTg0NCAyLjk2OTY4IDE0LjMwMzIgMS44NjA5NCAxMS40NDA5WiIvPgogICAgPHBhdGggY2xhc3M9ImpwLWljb24yIiBzdHJva2U9IiMzMzMzMzMiIHN0cm9rZS13aWR0aD0iMiIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoOS4zMTU5MiA5LjMyMDMxKSIgZD0iTTcuMzY4NDIgMEwwIDcuMzY0NzkiLz4KICAgIDxwYXRoIGNsYXNzPSJqcC1pY29uMiIgc3Ryb2tlPSIjMzMzMzMzIiBzdHJva2Utd2lkdGg9IjIiIHRyYW5zZm9ybT0idHJhbnNsYXRlKDkuMzE1OTIgMTYuNjgzNikgc2NhbGUoMSAtMSkiIGQ9Ik03LjM2ODQyIDBMMCA3LjM2NDc5Ii8+Cjwvc3ZnPgo=);
  --jp-icon-notebook: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8ZyBjbGFzcz0ianAtaWNvbi13YXJuMCBqcC1pY29uLXNlbGVjdGFibGUiIGZpbGw9IiNFRjZDMDAiPgogICAgPHBhdGggZD0iTTE4LjcgMy4zdjE1LjRIMy4zVjMuM2gxNS40bTEuNS0xLjVIMS44djE4LjNoMTguM2wuMS0xOC4zeiIvPgogICAgPHBhdGggZD0iTTE2LjUgMTYuNWwtNS40LTQuMy01LjYgNC4zdi0xMWgxMXoiLz4KICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-palette: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTE4IDEzVjIwSDRWNkg5LjAyQzkuMDcgNS4yOSA5LjI0IDQuNjIgOS41IDRINEMyLjkgNCAyIDQuOSAyIDZWMjBDMiAyMS4xIDIuOSAyMiA0IDIySDE4QzE5LjEgMjIgMjAgMjEuMSAyMCAyMFYxNUwxOCAxM1pNMTkuMyA4Ljg5QzE5Ljc0IDguMTkgMjAgNy4zOCAyMCA2LjVDMjAgNC4wMSAxNy45OSAyIDE1LjUgMkMxMy4wMSAyIDExIDQuMDEgMTEgNi41QzExIDguOTkgMTMuMDEgMTEgMTUuNDkgMTFDMTYuMzcgMTEgMTcuMTkgMTAuNzQgMTcuODggMTAuM0wyMSAxMy40MkwyMi40MiAxMkwxOS4zIDguODlaTTE1LjUgOUMxNC4xMiA5IDEzIDcuODggMTMgNi41QzEzIDUuMTIgMTQuMTIgNCAxNS41IDRDMTYuODggNCAxOCA1LjEyIDE4IDYuNUMxOCA3Ljg4IDE2Ljg4IDkgMTUuNSA5WiIvPgogICAgPHBhdGggZmlsbC1ydWxlPSJldmVub2RkIiBjbGlwLXJ1bGU9ImV2ZW5vZGQiIGQ9Ik00IDZIOS4wMTg5NEM5LjAwNjM5IDYuMTY1MDIgOSA2LjMzMTc2IDkgNi41QzkgOC44MTU3NyAxMC4yMTEgMTAuODQ4NyAxMi4wMzQzIDEySDlWMTRIMTZWMTIuOTgxMUMxNi41NzAzIDEyLjkzNzcgMTcuMTIgMTIuODIwNyAxNy42Mzk2IDEyLjYzOTZMMTggMTNWMjBINFY2Wk04IDhINlYxMEg4VjhaTTYgMTJIOFYxNEg2VjEyWk04IDE2SDZWMThIOFYxNlpNOSAxNkgxNlYxOEg5VjE2WiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-paste: url(data:image/svg+xml;base64,PHN2ZyBoZWlnaHQ9IjI0IiB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICAgIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICAgICAgPHBhdGggZD0iTTE5IDJoLTQuMThDMTQuNC44NCAxMy4zIDAgMTIgMGMtMS4zIDAtMi40Ljg0LTIuODIgMkg1Yy0xLjEgMC0yIC45LTIgMnYxNmMwIDEuMS45IDIgMiAyaDE0YzEuMSAwIDItLjkgMi0yVjRjMC0xLjEtLjktMi0yLTJ6bS03IDBjLjU1IDAgMSAuNDUgMSAxcy0uNDUgMS0xIDEtMS0uNDUtMS0xIC40NS0xIDEtMXptNyAxOEg1VjRoMnYzaDEwVjRoMnYxNnoiLz4KICAgIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-python: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8ZyBjbGFzcz0ianAtaWNvbi1icmFuZDAganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjMEQ0N0ExIj4KICAgIDxwYXRoIGQ9Ik0xMS4xIDYuOVY1LjhINi45YzAtLjUgMC0xLjMuMi0xLjYuNC0uNy44LTEuMSAxLjctMS40IDEuNy0uMyAyLjUtLjMgMy45LS4xIDEgLjEgMS45LjkgMS45IDEuOXY0LjJjMCAuNS0uOSAxLjYtMiAxLjZIOC44Yy0xLjUgMC0yLjQgMS40LTIuNCAyLjh2Mi4ySDQuN0MzLjUgMTUuMSAzIDE0IDMgMTMuMVY5Yy0uMS0xIC42LTIgMS44LTIgMS41LS4xIDYuMy0uMSA2LjMtLjF6Ii8+CiAgICA8cGF0aCBkPSJNMTAuOSAxNS4xdjEuMWg0LjJjMCAuNSAwIDEuMy0uMiAxLjYtLjQuNy0uOCAxLjEtMS43IDEuNC0xLjcuMy0yLjUuMy0zLjkuMS0xLS4xLTEuOS0uOS0xLjktMS45di00LjJjMC0uNS45LTEuNiAyLTEuNmgzLjhjMS41IDAgMi40LTEuNCAyLjQtMi44VjYuNmgxLjdDMTguNSA2LjkgMTkgOCAxOSA4LjlWMTNjMCAxLS43IDIuMS0xLjkgMi4xaC02LjJ6Ii8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-r-kernel: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8cGF0aCBjbGFzcz0ianAtaWNvbi1jb250cmFzdDMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjMjE5NkYzIiBkPSJNNC40IDIuNWMxLjItLjEgMi45LS4zIDQuOS0uMyAyLjUgMCA0LjEuNCA1LjIgMS4zIDEgLjcgMS41IDEuOSAxLjUgMy41IDAgMi0xLjQgMy41LTIuOSA0LjEgMS4yLjQgMS43IDEuNiAyLjIgMyAuNiAxLjkgMSAzLjkgMS4zIDQuNmgtMy44Yy0uMy0uNC0uOC0xLjctMS4yLTMuN3MtMS4yLTIuNi0yLjYtMi42aC0uOXY2LjRINC40VjIuNXptMy43IDYuOWgxLjRjMS45IDAgMi45LS45IDIuOS0yLjNzLTEtMi4zLTIuOC0yLjNjLS43IDAtMS4zIDAtMS42LjJ2NC41aC4xdi0uMXoiLz4KPC9zdmc+Cg==);
  --jp-icon-react: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMTUwIDE1MCA1NDEuOSAyOTUuMyI+CiAgPGcgY2xhc3M9ImpwLWljb24tYnJhbmQyIGpwLWljb24tc2VsZWN0YWJsZSIgZmlsbD0iIzYxREFGQiI+CiAgICA8cGF0aCBkPSJNNjY2LjMgMjk2LjVjMC0zMi41LTQwLjctNjMuMy0xMDMuMS04Mi40IDE0LjQtNjMuNiA4LTExNC4yLTIwLjItMTMwLjQtNi41LTMuOC0xNC4xLTUuNi0yMi40LTUuNnYyMi4zYzQuNiAwIDguMy45IDExLjQgMi42IDEzLjYgNy44IDE5LjUgMzcuNSAxNC45IDc1LjctMS4xIDkuNC0yLjkgMTkuMy01LjEgMjkuNC0xOS42LTQuOC00MS04LjUtNjMuNS0xMC45LTEzLjUtMTguNS0yNy41LTM1LjMtNDEuNi01MCAzMi42LTMwLjMgNjMuMi00Ni45IDg0LTQ2LjlWNzhjLTI3LjUgMC02My41IDE5LjYtOTkuOSA1My42LTM2LjQtMzMuOC03Mi40LTUzLjItOTkuOS01My4ydjIyLjNjMjAuNyAwIDUxLjQgMTYuNSA4NCA0Ni42LTE0IDE0LjctMjggMzEuNC00MS4zIDQ5LjktMjIuNiAyLjQtNDQgNi4xLTYzLjYgMTEtMi4zLTEwLTQtMTkuNy01LjItMjktNC43LTM4LjIgMS4xLTY3LjkgMTQuNi03NS44IDMtMS44IDYuOS0yLjYgMTEuNS0yLjZWNzguNWMtOC40IDAtMTYgMS44LTIyLjYgNS42LTI4LjEgMTYuMi0zNC40IDY2LjctMTkuOSAxMzAuMS02Mi4yIDE5LjItMTAyLjcgNDkuOS0xMDIuNyA4Mi4zIDAgMzIuNSA0MC43IDYzLjMgMTAzLjEgODIuNC0xNC40IDYzLjYtOCAxMTQuMiAyMC4yIDEzMC40IDYuNSAzLjggMTQuMSA1LjYgMjIuNSA1LjYgMjcuNSAwIDYzLjUtMTkuNiA5OS45LTUzLjYgMzYuNCAzMy44IDcyLjQgNTMuMiA5OS45IDUzLjIgOC40IDAgMTYtMS44IDIyLjYtNS42IDI4LjEtMTYuMiAzNC40LTY2LjcgMTkuOS0xMzAuMSA2Mi0xOS4xIDEwMi41LTQ5LjkgMTAyLjUtODIuM3ptLTEzMC4yLTY2LjdjLTMuNyAxMi45LTguMyAyNi4yLTEzLjUgMzkuNS00LjEtOC04LjQtMTYtMTMuMS0yNC00LjYtOC05LjUtMTUuOC0xNC40LTIzLjQgMTQuMiAyLjEgMjcuOSA0LjcgNDEgNy45em0tNDUuOCAxMDYuNWMtNy44IDEzLjUtMTUuOCAyNi4zLTI0LjEgMzguMi0xNC45IDEuMy0zMCAyLTQ1LjIgMi0xNS4xIDAtMzAuMi0uNy00NS0xLjktOC4zLTExLjktMTYuNC0yNC42LTI0LjItMzgtNy42LTEzLjEtMTQuNS0yNi40LTIwLjgtMzkuOCA2LjItMTMuNCAxMy4yLTI2LjggMjAuNy0zOS45IDcuOC0xMy41IDE1LjgtMjYuMyAyNC4xLTM4LjIgMTQuOS0xLjMgMzAtMiA0NS4yLTIgMTUuMSAwIDMwLjIuNyA0NSAxLjkgOC4zIDExLjkgMTYuNCAyNC42IDI0LjIgMzggNy42IDEzLjEgMTQuNSAyNi40IDIwLjggMzkuOC02LjMgMTMuNC0xMy4yIDI2LjgtMjAuNyAzOS45em0zMi4zLTEzYzUuNCAxMy40IDEwIDI2LjggMTMuOCAzOS44LTEzLjEgMy4yLTI2LjkgNS45LTQxLjIgOCA0LjktNy43IDkuOC0xNS42IDE0LjQtMjMuNyA0LjYtOCA4LjktMTYuMSAxMy0yNC4xek00MjEuMiA0MzBjLTkuMy05LjYtMTguNi0yMC4zLTI3LjgtMzIgOSAuNCAxOC4yLjcgMjcuNS43IDkuNCAwIDE4LjctLjIgMjcuOC0uNy05IDExLjctMTguMyAyMi40LTI3LjUgMzJ6bS03NC40LTU4LjljLTE0LjItMi4xLTI3LjktNC43LTQxLTcuOSAzLjctMTIuOSA4LjMtMjYuMiAxMy41LTM5LjUgNC4xIDggOC40IDE2IDEzLjEgMjQgNC43IDggOS41IDE1LjggMTQuNCAyMy40ek00MjAuNyAxNjNjOS4zIDkuNiAxOC42IDIwLjMgMjcuOCAzMi05LS40LTE4LjItLjctMjcuNS0uNy05LjQgMC0xOC43LjItMjcuOC43IDktMTEuNyAxOC4zLTIyLjQgMjcuNS0zMnptLTc0IDU4LjljLTQuOSA3LjctOS44IDE1LjYtMTQuNCAyMy43LTQuNiA4LTguOSAxNi0xMyAyNC01LjQtMTMuNC0xMC0yNi44LTEzLjgtMzkuOCAxMy4xLTMuMSAyNi45LTUuOCA0MS4yLTcuOXptLTkwLjUgMTI1LjJjLTM1LjQtMTUuMS01OC4zLTM0LjktNTguMy01MC42IDAtMTUuNyAyMi45LTM1LjYgNTguMy01MC42IDguNi0zLjcgMTgtNyAyNy43LTEwLjEgNS43IDE5LjYgMTMuMiA0MCAyMi41IDYwLjktOS4yIDIwLjgtMTYuNiA0MS4xLTIyLjIgNjAuNi05LjktMy4xLTE5LjMtNi41LTI4LTEwLjJ6TTMxMCA0OTBjLTEzLjYtNy44LTE5LjUtMzcuNS0xNC45LTc1LjcgMS4xLTkuNCAyLjktMTkuMyA1LjEtMjkuNCAxOS42IDQuOCA0MSA4LjUgNjMuNSAxMC45IDEzLjUgMTguNSAyNy41IDM1LjMgNDEuNiA1MC0zMi42IDMwLjMtNjMuMiA0Ni45LTg0IDQ2LjktNC41LS4xLTguMy0xLTExLjMtMi43em0yMzcuMi03Ni4yYzQuNyAzOC4yLTEuMSA2Ny45LTE0LjYgNzUuOC0zIDEuOC02LjkgMi42LTExLjUgMi42LTIwLjcgMC01MS40LTE2LjUtODQtNDYuNiAxNC0xNC43IDI4LTMxLjQgNDEuMy00OS45IDIyLjYtMi40IDQ0LTYuMSA2My42LTExIDIuMyAxMC4xIDQuMSAxOS44IDUuMiAyOS4xem0zOC41LTY2LjdjLTguNiAzLjctMTggNy0yNy43IDEwLjEtNS43LTE5LjYtMTMuMi00MC0yMi41LTYwLjkgOS4yLTIwLjggMTYuNi00MS4xIDIyLjItNjAuNiA5LjkgMy4xIDE5LjMgNi41IDI4LjEgMTAuMiAzNS40IDE1LjEgNTguMyAzNC45IDU4LjMgNTAuNi0uMSAxNS43LTIzIDM1LjYtNTguNCA1MC42ek0zMjAuOCA3OC40eiIvPgogICAgPGNpcmNsZSBjeD0iNDIwLjkiIGN5PSIyOTYuNSIgcj0iNDUuNyIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-refresh: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDE4IDE4Ij4KICAgIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICAgICAgPHBhdGggZD0iTTkgMTMuNWMtMi40OSAwLTQuNS0yLjAxLTQuNS00LjVTNi41MSA0LjUgOSA0LjVjMS4yNCAwIDIuMzYuNTIgMy4xNyAxLjMzTDEwIDhoNVYzbC0xLjc2IDEuNzZDMTIuMTUgMy42OCAxMC42NiAzIDkgMyA1LjY5IDMgMy4wMSA1LjY5IDMuMDEgOVM1LjY5IDE1IDkgMTVjMi45NyAwIDUuNDMtMi4xNiA1LjktNWgtMS41MmMtLjQ2IDItMi4yNCAzLjUtNC4zOCAzLjV6Ii8+CiAgICA8L2c+Cjwvc3ZnPgo=);
  --jp-icon-regex: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIwIDIwIj4KICA8ZyBjbGFzcz0ianAtaWNvbjIiIGZpbGw9IiM0MTQxNDEiPgogICAgPHJlY3QgeD0iMiIgeT0iMiIgd2lkdGg9IjE2IiBoZWlnaHQ9IjE2Ii8+CiAgPC9nPgoKICA8ZyBjbGFzcz0ianAtaWNvbi1hY2NlbnQyIiBmaWxsPSIjRkZGIj4KICAgIDxjaXJjbGUgY2xhc3M9InN0MiIgY3g9IjUuNSIgY3k9IjE0LjUiIHI9IjEuNSIvPgogICAgPHJlY3QgeD0iMTIiIHk9IjQiIGNsYXNzPSJzdDIiIHdpZHRoPSIxIiBoZWlnaHQ9IjgiLz4KICAgIDxyZWN0IHg9IjguNSIgeT0iNy41IiB0cmFuc2Zvcm09Im1hdHJpeCgwLjg2NiAtMC41IDAuNSAwLjg2NiAtMi4zMjU1IDcuMzIxOSkiIGNsYXNzPSJzdDIiIHdpZHRoPSI4IiBoZWlnaHQ9IjEiLz4KICAgIDxyZWN0IHg9IjEyIiB5PSI0IiB0cmFuc2Zvcm09Im1hdHJpeCgwLjUgLTAuODY2IDAuODY2IDAuNSAtMC42Nzc5IDE0LjgyNTIpIiBjbGFzcz0ic3QyIiB3aWR0aD0iMSIgaGVpZ2h0PSI4Ii8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-run: url(data:image/svg+xml;base64,PHN2ZyBoZWlnaHQ9IjI0IiB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICAgIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICAgICAgPHBhdGggZD0iTTggNXYxNGwxMS03eiIvPgogICAgPC9nPgo8L3N2Zz4K);
  --jp-icon-running: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDUxMiA1MTIiPgogIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICA8cGF0aCBkPSJNMjU2IDhDMTE5IDggOCAxMTkgOCAyNTZzMTExIDI0OCAyNDggMjQ4IDI0OC0xMTEgMjQ4LTI0OFMzOTMgOCAyNTYgOHptOTYgMzI4YzAgOC44LTcuMiAxNi0xNiAxNkgxNzZjLTguOCAwLTE2LTcuMi0xNi0xNlYxNzZjMC04LjggNy4yLTE2IDE2LTE2aDE2MGM4LjggMCAxNiA3LjIgMTYgMTZ2MTYweiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-save: url(data:image/svg+xml;base64,PHN2ZyBoZWlnaHQ9IjI0IiB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICAgIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICAgICAgPHBhdGggZD0iTTE3IDNINWMtMS4xMSAwLTIgLjktMiAydjE0YzAgMS4xLjg5IDIgMiAyaDE0YzEuMSAwIDItLjkgMi0yVjdsLTQtNHptLTUgMTZjLTEuNjYgMC0zLTEuMzQtMy0zczEuMzQtMyAzLTMgMyAxLjM0IDMgMy0xLjM0IDMtMyAzem0zLTEwSDVWNWgxMHY0eiIvPgogICAgPC9nPgo8L3N2Zz4K);
  --jp-icon-search: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMTggMTgiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTEyLjEsMTAuOWgtMC43bC0wLjItMC4yYzAuOC0wLjksMS4zLTIuMiwxLjMtMy41YzAtMy0yLjQtNS40LTUuNC01LjRTMS44LDQuMiwxLjgsNy4xczIuNCw1LjQsNS40LDUuNCBjMS4zLDAsMi41LTAuNSwzLjUtMS4zbDAuMiwwLjJ2MC43bDQuMSw0LjFsMS4yLTEuMkwxMi4xLDEwLjl6IE03LjEsMTAuOWMtMi4xLDAtMy43LTEuNy0zLjctMy43czEuNy0zLjcsMy43LTMuN3MzLjcsMS43LDMuNywzLjcgUzkuMiwxMC45LDcuMSwxMC45eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-settings: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMTkuNDMgMTIuOThjLjA0LS4zMi4wNy0uNjQuMDctLjk4cy0uMDMtLjY2LS4wNy0uOThsMi4xMS0xLjY1Yy4xOS0uMTUuMjQtLjQyLjEyLS42NGwtMi0zLjQ2Yy0uMTItLjIyLS4zOS0uMy0uNjEtLjIybC0yLjQ5IDFjLS41Mi0uNC0xLjA4LS43My0xLjY5LS45OGwtLjM4LTIuNjVBLjQ4OC40ODggMCAwMDE0IDJoLTRjLS4yNSAwLS40Ni4xOC0uNDkuNDJsLS4zOCAyLjY1Yy0uNjEuMjUtMS4xNy41OS0xLjY5Ljk4bC0yLjQ5LTFjLS4yMy0uMDktLjQ5IDAtLjYxLjIybC0yIDMuNDZjLS4xMy4yMi0uMDcuNDkuMTIuNjRsMi4xMSAxLjY1Yy0uMDQuMzItLjA3LjY1LS4wNy45OHMuMDMuNjYuMDcuOThsLTIuMTEgMS42NWMtLjE5LjE1LS4yNC40Mi0uMTIuNjRsMiAzLjQ2Yy4xMi4yMi4zOS4zLjYxLjIybDIuNDktMWMuNTIuNCAxLjA4LjczIDEuNjkuOThsLjM4IDIuNjVjLjAzLjI0LjI0LjQyLjQ5LjQyaDRjLjI1IDAgLjQ2LS4xOC40OS0uNDJsLjM4LTIuNjVjLjYxLS4yNSAxLjE3LS41OSAxLjY5LS45OGwyLjQ5IDFjLjIzLjA5LjQ5IDAgLjYxLS4yMmwyLTMuNDZjLjEyLS4yMi4wNy0uNDktLjEyLS42NGwtMi4xMS0xLjY1ek0xMiAxNS41Yy0xLjkzIDAtMy41LTEuNTctMy41LTMuNXMxLjU3LTMuNSAzLjUtMy41IDMuNSAxLjU3IDMuNSAzLjUtMS41NyAzLjUtMy41IDMuNXoiLz4KPC9zdmc+Cg==);
  --jp-icon-spreadsheet: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8cGF0aCBjbGFzcz0ianAtaWNvbi1jb250cmFzdDEganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNENBRjUwIiBkPSJNMi4yIDIuMnYxNy42aDE3LjZWMi4ySDIuMnptMTUuNCA3LjdoLTUuNVY0LjRoNS41djUuNXpNOS45IDQuNHY1LjVINC40VjQuNGg1LjV6bS01LjUgNy43aDUuNXY1LjVINC40di01LjV6bTcuNyA1LjV2LTUuNWg1LjV2NS41aC01LjV6Ii8+Cjwvc3ZnPgo=);
  --jp-icon-stop: url(data:image/svg+xml;base64,PHN2ZyBoZWlnaHQ9IjI0IiB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICAgIDxnIGNsYXNzPSJqcC1pY29uMyIgZmlsbD0iIzYxNjE2MSI+CiAgICAgICAgPHBhdGggZD0iTTAgMGgyNHYyNEgweiIgZmlsbD0ibm9uZSIvPgogICAgICAgIDxwYXRoIGQ9Ik02IDZoMTJ2MTJINnoiLz4KICAgIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-tab: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTIxIDNIM2MtMS4xIDAtMiAuOS0yIDJ2MTRjMCAxLjEuOSAyIDIgMmgxOGMxLjEgMCAyLS45IDItMlY1YzAtMS4xLS45LTItMi0yem0wIDE2SDNWNWgxMHY0aDh2MTB6Ii8+CiAgPC9nPgo8L3N2Zz4K);
  --jp-icon-terminal: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0IiA+CiAgICA8cmVjdCBjbGFzcz0ianAtaWNvbjIganAtaWNvbi1zZWxlY3RhYmxlIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHRyYW5zZm9ybT0idHJhbnNsYXRlKDIgMikiIGZpbGw9IiMzMzMzMzMiLz4KICAgIDxwYXRoIGNsYXNzPSJqcC1pY29uLWFjY2VudDIganAtaWNvbi1zZWxlY3RhYmxlLWludmVyc2UiIGQ9Ik01LjA1NjY0IDguNzYxNzJDNS4wNTY2NCA4LjU5NzY2IDUuMDMxMjUgOC40NTMxMiA0Ljk4MDQ3IDguMzI4MTJDNC45MzM1OSA4LjE5OTIyIDQuODU1NDcgOC4wODIwMyA0Ljc0NjA5IDcuOTc2NTZDNC42NDA2MiA3Ljg3MTA5IDQuNSA3Ljc3NTM5IDQuMzI0MjIgNy42ODk0NUM0LjE1MjM0IDcuNTk5NjEgMy45NDMzNiA3LjUxMTcyIDMuNjk3MjcgNy40MjU3OEMzLjMwMjczIDcuMjg1MTYgMi45NDMzNiA3LjEzNjcyIDIuNjE5MTQgNi45ODA0N0MyLjI5NDkyIDYuODI0MjIgMi4wMTc1OCA2LjY0MjU4IDEuNzg3MTEgNi40MzU1NUMxLjU2MDU1IDYuMjI4NTIgMS4zODQ3NyA1Ljk4ODI4IDEuMjU5NzcgNS43MTQ4NEMxLjEzNDc3IDUuNDM3NSAxLjA3MjI3IDUuMTA5MzggMS4wNzIyNyA0LjczMDQ3QzEuMDcyMjcgNC4zOTg0NCAxLjEyODkxIDQuMDk1NyAxLjI0MjE5IDMuODIyMjdDMS4zNTU0NyAzLjU0NDkyIDEuNTE1NjIgMy4zMDQ2OSAxLjcyMjY2IDMuMTAxNTZDMS45Mjk2OSAyLjg5ODQ0IDIuMTc5NjkgMi43MzQzNyAyLjQ3MjY2IDIuNjA5MzhDMi43NjU2MiAyLjQ4NDM4IDMuMDkxOCAyLjQwNDMgMy40NTExNyAyLjM2OTE0VjEuMTA5MzhINC4zODg2N1YyLjM4MDg2QzQuNzQwMjMgMi40Mjc3MyA1LjA1NjY0IDIuNTIzNDQgNS4zMzc4OSAyLjY2Nzk3QzUuNjE5MTQgMi44MTI1IDUuODU3NDIgMy4wMDE5NSA2LjA1MjczIDMuMjM2MzNDNi4yNTE5NSAzLjQ2NjggNi40MDQzIDMuNzQwMjMgNi41MDk3NyA0LjA1NjY0QzYuNjE5MTQgNC4zNjkxNCA2LjY3MzgzIDQuNzIwNyA2LjY3MzgzIDUuMTExMzNINS4wNDQ5MkM1LjA0NDkyIDQuNjM4NjcgNC45Mzc1IDQuMjgxMjUgNC43MjI2NiA0LjAzOTA2QzQuNTA3ODEgMy43OTI5NyA0LjIxNjggMy42Njk5MiAzLjg0OTYxIDMuNjY5OTJDMy42NTAzOSAzLjY2OTkyIDMuNDc2NTYgMy42OTcyNyAzLjMyODEyIDMuNzUxOTVDMy4xODM1OSAzLjgwMjczIDMuMDY0NDUgMy44NzY5NSAyLjk3MDcgMy45NzQ2MUMyLjg3Njk1IDQuMDY4MzYgMi44MDY2NCA0LjE3OTY5IDIuNzU5NzcgNC4zMDg1OUMyLjcxNjggNC40Mzc1IDIuNjk1MzEgNC41NzgxMiAyLjY5NTMxIDQuNzMwNDdDMi42OTUzMSA0Ljg4MjgxIDIuNzE2OCA1LjAxOTUzIDIuNzU5NzcgNS4xNDA2MkMyLjgwNjY0IDUuMjU3ODEgMi44ODI4MSA1LjM2NzE5IDIuOTg4MjggNS40Njg3NUMzLjA5NzY2IDUuNTcwMzEgMy4yNDAyMyA1LjY2Nzk3IDMuNDE2MDIgNS43NjE3MkMzLjU5MTggNS44NTE1NiAzLjgxMDU1IDUuOTQzMzYgNC4wNzIyNyA2LjAzNzExQzQuNDY2OCA2LjE4NTU1IDQuODI0MjIgNi4zMzk4NCA1LjE0NDUzIDYuNUM1LjQ2NDg0IDYuNjU2MjUgNS43MzgyOCA2LjgzOTg0IDUuOTY0ODQgNy4wNTA3OEM2LjE5NTMxIDcuMjU3ODEgNi4zNzEwOSA3LjUgNi40OTIxOSA3Ljc3NzM0QzYuNjE3MTkgOC4wNTA3OCA2LjY3OTY5IDguMzc1IDYuNjc5NjkgOC43NUM2LjY3OTY5IDkuMDkzNzUgNi42MjMwNSA5LjQwNDMgNi41MDk3NyA5LjY4MTY0QzYuMzk2NDggOS45NTUwOCA2LjIzNDM4IDEwLjE5MTQgNi4wMjM0NCAxMC4zOTA2QzUuODEyNSAxMC41ODk4IDUuNTU4NTkgMTAuNzUgNS4yNjE3MiAxMC44NzExQzQuOTY0ODQgMTAuOTg4MyA0LjYzMjgxIDExLjA2NDUgNC4yNjU2MiAxMS4wOTk2VjEyLjI0OEgzLjMzMzk4VjExLjA5OTZDMy4wMDE5NSAxMS4wNjg0IDIuNjc5NjkgMTAuOTk2MSAyLjM2NzE5IDEwLjg4MjhDMi4wNTQ2OSAxMC43NjU2IDEuNzc3MzQgMTAuNTk3NyAxLjUzNTE2IDEwLjM3ODlDMS4yOTY4OCAxMC4xNjAyIDEuMTA1NDcgOS44ODQ3NyAwLjk2MDkzOCA5LjU1MjczQzAuODE2NDA2IDkuMjE2OCAwLjc0NDE0MSA4LjgxNDQ1IDAuNzQ0MTQxIDguMzQ1N0gyLjM3ODkxQzIuMzc4OTEgOC42MjY5NSAyLjQxOTkyIDguODYzMjggMi41MDE5NSA5LjA1NDY5QzIuNTgzOTggOS4yNDIxOSAyLjY4OTQ1IDkuMzkyNTggMi44MTgzNiA5LjUwNTg2QzIuOTUxMTcgOS42MTUyMyAzLjEwMTU2IDkuNjkzMzYgMy4yNjk1MyA5Ljc0MDIzQzMuNDM3NSA5Ljc4NzExIDMuNjA5MzggOS44MTA1NSAzLjc4NTE2IDkuODEwNTVDNC4yMDMxMiA5LjgxMDU1IDQuNTE5NTMgOS43MTI4OSA0LjczNDM4IDkuNTE3NThDNC45NDkyMiA5LjMyMjI3IDUuMDU2NjQgOS4wNzAzMSA1LjA1NjY0IDguNzYxNzJaTTEzLjQxOCAxMi4yNzE1SDguMDc0MjJWMTFIMTMuNDE4VjEyLjI3MTVaIiB0cmFuc2Zvcm09InRyYW5zbGF0ZSgzLjk1MjY0IDYpIiBmaWxsPSJ3aGl0ZSIvPgo8L3N2Zz4K);
  --jp-icon-text-editor: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI0Ij4KICA8cGF0aCBjbGFzcz0ianAtaWNvbjMganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjNjE2MTYxIiBkPSJNMTUgMTVIM3YyaDEydi0yem0wLThIM3YyaDEyVjd6TTMgMTNoMTh2LTJIM3Yyem0wIDhoMTh2LTJIM3Yyek0zIDN2MmgxOFYzSDN6Ii8+Cjwvc3ZnPgo=);
  --jp-icon-trusted: url(data:image/svg+xml;base64,PHN2ZyBmaWxsPSJub25lIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDI0IDI1Ij4KICAgIDxwYXRoIGNsYXNzPSJqcC1pY29uMiIgc3Ryb2tlPSIjMzMzMzMzIiBzdHJva2Utd2lkdGg9IjIiIHRyYW5zZm9ybT0idHJhbnNsYXRlKDIgMykiIGQ9Ik0xLjg2MDk0IDExLjQ0MDlDMC44MjY0NDggOC43NzAyNyAwLjg2Mzc3OSA2LjA1NzY0IDEuMjQ5MDcgNC4xOTkzMkMyLjQ4MjA2IDMuOTMzNDcgNC4wODA2OCAzLjQwMzQ3IDUuNjAxMDIgMi44NDQ5QzcuMjM1NDkgMi4yNDQ0IDguODU2NjYgMS41ODE1IDkuOTg3NiAxLjA5NTM5QzExLjA1OTcgMS41ODM0MSAxMi42MDk0IDIuMjQ0NCAxNC4yMTggMi44NDMzOUMxNS43NTAzIDMuNDEzOTQgMTcuMzk5NSAzLjk1MjU4IDE4Ljc1MzkgNC4yMTM4NUMxOS4xMzY0IDYuMDcxNzcgMTkuMTcwOSA4Ljc3NzIyIDE4LjEzOSAxMS40NDA5QzE3LjAzMDMgMTQuMzAzMiAxNC42NjY4IDE3LjE4NDQgOS45OTk5OSAxOC45MzU0QzUuMzMzMiAxNy4xODQ0IDIuOTY5NjggMTQuMzAzMiAxLjg2MDk0IDExLjQ0MDlaIi8+CiAgICA8cGF0aCBjbGFzcz0ianAtaWNvbjIiIGZpbGw9IiMzMzMzMzMiIHN0cm9rZT0iIzMzMzMzMyIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoOCA5Ljg2NzE5KSIgZD0iTTIuODYwMTUgNC44NjUzNUwwLjcyNjU0OSAyLjk5OTU5TDAgMy42MzA0NUwyLjg2MDE1IDYuMTMxNTdMOCAwLjYzMDg3Mkw3LjI3ODU3IDBMMi44NjAxNSA0Ljg2NTM1WiIvPgo8L3N2Zz4K);
  --jp-icon-undo: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMjQgMjQiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTEyLjUgOGMtMi42NSAwLTUuMDUuOTktNi45IDIuNkwyIDd2OWg5bC0zLjYyLTMuNjJjMS4zOS0xLjE2IDMuMTYtMS44OCA1LjEyLTEuODggMy41NCAwIDYuNTUgMi4zMSA3LjYgNS41bDIuMzctLjc4QzIxLjA4IDExLjAzIDE3LjE1IDggMTIuNSA4eiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-vega: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8ZyBjbGFzcz0ianAtaWNvbjEganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjMjEyMTIxIj4KICAgIDxwYXRoIGQ9Ik0xMC42IDUuNGwyLjItMy4ySDIuMnY3LjNsNC02LjZ6Ii8+CiAgICA8cGF0aCBkPSJNMTUuOCAyLjJsLTQuNCA2LjZMNyA2LjNsLTQuOCA4djUuNWgxNy42VjIuMmgtNHptLTcgMTUuNEg1LjV2LTQuNGgzLjN2NC40em00LjQgMEg5LjhWOS44aDMuNHY3Ljh6bTQuNCAwaC0zLjRWNi41aDMuNHYxMS4xeiIvPgogIDwvZz4KPC9zdmc+Cg==);
  --jp-icon-yaml: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgdmlld0JveD0iMCAwIDIyIDIyIj4KICA8ZyBjbGFzcz0ianAtaWNvbi1jb250cmFzdDIganAtaWNvbi1zZWxlY3RhYmxlIiBmaWxsPSIjRDgxQjYwIj4KICAgIDxwYXRoIGQ9Ik03LjIgMTguNnYtNS40TDMgNS42aDMuM2wxLjQgMy4xYy4zLjkuNiAxLjYgMSAyLjUuMy0uOC42LTEuNiAxLTIuNWwxLjQtMy4xaDMuNGwtNC40IDcuNnY1LjVsLTIuOS0uMXoiLz4KICAgIDxjaXJjbGUgY2xhc3M9InN0MCIgY3g9IjE3LjYiIGN5PSIxNi41IiByPSIyLjEiLz4KICAgIDxjaXJjbGUgY2xhc3M9InN0MCIgY3g9IjE3LjYiIGN5PSIxMSIgcj0iMi4xIi8+CiAgPC9nPgo8L3N2Zz4K);
}

/* Icon CSS class declarations */

.jp-AddIcon {
  background-image: var(--jp-icon-add);
}
.jp-BugIcon {
  background-image: var(--jp-icon-bug);
}
.jp-BuildIcon {
  background-image: var(--jp-icon-build);
}
.jp-CaretDownEmptyIcon {
  background-image: var(--jp-icon-caret-down-empty);
}
.jp-CaretDownEmptyThinIcon {
  background-image: var(--jp-icon-caret-down-empty-thin);
}
.jp-CaretDownIcon {
  background-image: var(--jp-icon-caret-down);
}
.jp-CaretLeftIcon {
  background-image: var(--jp-icon-caret-left);
}
.jp-CaretRightIcon {
  background-image: var(--jp-icon-caret-right);
}
.jp-CaretUpEmptyThinIcon {
  background-image: var(--jp-icon-caret-up-empty-thin);
}
.jp-CaretUpIcon {
  background-image: var(--jp-icon-caret-up);
}
.jp-CaseSensitiveIcon {
  background-image: var(--jp-icon-case-sensitive);
}
.jp-CheckIcon {
  background-image: var(--jp-icon-check);
}
.jp-CircleEmptyIcon {
  background-image: var(--jp-icon-circle-empty);
}
.jp-CircleIcon {
  background-image: var(--jp-icon-circle);
}
.jp-ClearIcon {
  background-image: var(--jp-icon-clear);
}
.jp-CloseIcon {
  background-image: var(--jp-icon-close);
}
.jp-ConsoleIcon {
  background-image: var(--jp-icon-console);
}
.jp-CopyIcon {
  background-image: var(--jp-icon-copy);
}
.jp-CutIcon {
  background-image: var(--jp-icon-cut);
}
.jp-DownloadIcon {
  background-image: var(--jp-icon-download);
}
.jp-EditIcon {
  background-image: var(--jp-icon-edit);
}
.jp-EllipsesIcon {
  background-image: var(--jp-icon-ellipses);
}
.jp-ExtensionIcon {
  background-image: var(--jp-icon-extension);
}
.jp-FastForwardIcon {
  background-image: var(--jp-icon-fast-forward);
}
.jp-FileIcon {
  background-image: var(--jp-icon-file);
}
.jp-FileUploadIcon {
  background-image: var(--jp-icon-file-upload);
}
.jp-FilterListIcon {
  background-image: var(--jp-icon-filter-list);
}
.jp-FolderIcon {
  background-image: var(--jp-icon-folder);
}
.jp-Html5Icon {
  background-image: var(--jp-icon-html5);
}
.jp-ImageIcon {
  background-image: var(--jp-icon-image);
}
.jp-InspectorIcon {
  background-image: var(--jp-icon-inspector);
}
.jp-JsonIcon {
  background-image: var(--jp-icon-json);
}
.jp-JupyterFaviconIcon {
  background-image: var(--jp-icon-jupyter-favicon);
}
.jp-JupyterIcon {
  background-image: var(--jp-icon-jupyter);
}
.jp-JupyterlabWordmarkIcon {
  background-image: var(--jp-icon-jupyterlab-wordmark);
}
.jp-KernelIcon {
  background-image: var(--jp-icon-kernel);
}
.jp-KeyboardIcon {
  background-image: var(--jp-icon-keyboard);
}
.jp-LauncherIcon {
  background-image: var(--jp-icon-launcher);
}
.jp-LineFormIcon {
  background-image: var(--jp-icon-line-form);
}
.jp-LinkIcon {
  background-image: var(--jp-icon-link);
}
.jp-ListIcon {
  background-image: var(--jp-icon-list);
}
.jp-ListingsInfoIcon {
  background-image: var(--jp-icon-listings-info);
}
.jp-MarkdownIcon {
  background-image: var(--jp-icon-markdown);
}
.jp-NewFolderIcon {
  background-image: var(--jp-icon-new-folder);
}
.jp-NotTrustedIcon {
  background-image: var(--jp-icon-not-trusted);
}
.jp-NotebookIcon {
  background-image: var(--jp-icon-notebook);
}
.jp-PaletteIcon {
  background-image: var(--jp-icon-palette);
}
.jp-PasteIcon {
  background-image: var(--jp-icon-paste);
}
.jp-PythonIcon {
  background-image: var(--jp-icon-python);
}
.jp-RKernelIcon {
  background-image: var(--jp-icon-r-kernel);
}
.jp-ReactIcon {
  background-image: var(--jp-icon-react);
}
.jp-RefreshIcon {
  background-image: var(--jp-icon-refresh);
}
.jp-RegexIcon {
  background-image: var(--jp-icon-regex);
}
.jp-RunIcon {
  background-image: var(--jp-icon-run);
}
.jp-RunningIcon {
  background-image: var(--jp-icon-running);
}
.jp-SaveIcon {
  background-image: var(--jp-icon-save);
}
.jp-SearchIcon {
  background-image: var(--jp-icon-search);
}
.jp-SettingsIcon {
  background-image: var(--jp-icon-settings);
}
.jp-SpreadsheetIcon {
  background-image: var(--jp-icon-spreadsheet);
}
.jp-StopIcon {
  background-image: var(--jp-icon-stop);
}
.jp-TabIcon {
  background-image: var(--jp-icon-tab);
}
.jp-TerminalIcon {
  background-image: var(--jp-icon-terminal);
}
.jp-TextEditorIcon {
  background-image: var(--jp-icon-text-editor);
}
.jp-TrustedIcon {
  background-image: var(--jp-icon-trusted);
}
.jp-UndoIcon {
  background-image: var(--jp-icon-undo);
}
.jp-VegaIcon {
  background-image: var(--jp-icon-vega);
}
.jp-YamlIcon {
  background-image: var(--jp-icon-yaml);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/**
 * (DEPRECATED) Support for consuming icons as CSS background images
 */

:root {
  --jp-icon-search-white: url(data:image/svg+xml;base64,PHN2ZyB2aWV3Qm94PSIwIDAgMTggMTgiIHdpZHRoPSIxNiIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8ZyBjbGFzcz0ianAtaWNvbjMiIGZpbGw9IiM2MTYxNjEiPgogICAgPHBhdGggZD0iTTEyLjEsMTAuOWgtMC43bC0wLjItMC4yYzAuOC0wLjksMS4zLTIuMiwxLjMtMy41YzAtMy0yLjQtNS40LTUuNC01LjRTMS44LDQuMiwxLjgsNy4xczIuNCw1LjQsNS40LDUuNCBjMS4zLDAsMi41LTAuNSwzLjUtMS4zbDAuMiwwLjJ2MC43bDQuMSw0LjFsMS4yLTEuMkwxMi4xLDEwLjl6IE03LjEsMTAuOWMtMi4xLDAtMy43LTEuNy0zLjctMy43czEuNy0zLjcsMy43LTMuN3MzLjcsMS43LDMuNywzLjcgUzkuMiwxMC45LDcuMSwxMC45eiIvPgogIDwvZz4KPC9zdmc+Cg==);
}

.jp-Icon,
.jp-MaterialIcon {
  background-position: center;
  background-repeat: no-repeat;
  background-size: 16px;
  min-width: 16px;
  min-height: 16px;
}

.jp-Icon-cover {
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
}

/**
 * (DEPRECATED) Support for specific CSS icon sizes
 */

.jp-Icon-16 {
  background-size: 16px;
  min-width: 16px;
  min-height: 16px;
}

.jp-Icon-18 {
  background-size: 18px;
  min-width: 18px;
  min-height: 18px;
}

.jp-Icon-20 {
  background-size: 20px;
  min-width: 20px;
  min-height: 20px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/**
 * Support for icons as inline SVG HTMLElements
 */

/* recolor the primary elements of an icon */
.jp-icon0[fill] {
  fill: var(--jp-inverse-layout-color0);
}
.jp-icon1[fill] {
  fill: var(--jp-inverse-layout-color1);
}
.jp-icon2[fill] {
  fill: var(--jp-inverse-layout-color2);
}
.jp-icon3[fill] {
  fill: var(--jp-inverse-layout-color3);
}
.jp-icon4[fill] {
  fill: var(--jp-inverse-layout-color4);
}

.jp-icon0[stroke] {
  stroke: var(--jp-inverse-layout-color0);
}
.jp-icon1[stroke] {
  stroke: var(--jp-inverse-layout-color1);
}
.jp-icon2[stroke] {
  stroke: var(--jp-inverse-layout-color2);
}
.jp-icon3[stroke] {
  stroke: var(--jp-inverse-layout-color3);
}
.jp-icon4[stroke] {
  stroke: var(--jp-inverse-layout-color4);
}
/* recolor the accent elements of an icon */
.jp-icon-accent0[fill] {
  fill: var(--jp-layout-color0);
}
.jp-icon-accent1[fill] {
  fill: var(--jp-layout-color1);
}
.jp-icon-accent2[fill] {
  fill: var(--jp-layout-color2);
}
.jp-icon-accent3[fill] {
  fill: var(--jp-layout-color3);
}
.jp-icon-accent4[fill] {
  fill: var(--jp-layout-color4);
}

.jp-icon-accent0[stroke] {
  stroke: var(--jp-layout-color0);
}
.jp-icon-accent1[stroke] {
  stroke: var(--jp-layout-color1);
}
.jp-icon-accent2[stroke] {
  stroke: var(--jp-layout-color2);
}
.jp-icon-accent3[stroke] {
  stroke: var(--jp-layout-color3);
}
.jp-icon-accent4[stroke] {
  stroke: var(--jp-layout-color4);
}
/* set the color of an icon to transparent */
.jp-icon-none[fill] {
  fill: none;
}

.jp-icon-none[stroke] {
  stroke: none;
}
/* brand icon colors. Same for light and dark */
.jp-icon-brand0[fill] {
  fill: var(--jp-brand-color0);
}
.jp-icon-brand1[fill] {
  fill: var(--jp-brand-color1);
}
.jp-icon-brand2[fill] {
  fill: var(--jp-brand-color2);
}
.jp-icon-brand3[fill] {
  fill: var(--jp-brand-color3);
}
.jp-icon-brand4[fill] {
  fill: var(--jp-brand-color4);
}

.jp-icon-brand0[stroke] {
  stroke: var(--jp-brand-color0);
}
.jp-icon-brand1[stroke] {
  stroke: var(--jp-brand-color1);
}
.jp-icon-brand2[stroke] {
  stroke: var(--jp-brand-color2);
}
.jp-icon-brand3[stroke] {
  stroke: var(--jp-brand-color3);
}
.jp-icon-brand4[stroke] {
  stroke: var(--jp-brand-color4);
}
/* warn icon colors. Same for light and dark */
.jp-icon-warn0[fill] {
  fill: var(--jp-warn-color0);
}
.jp-icon-warn1[fill] {
  fill: var(--jp-warn-color1);
}
.jp-icon-warn2[fill] {
  fill: var(--jp-warn-color2);
}
.jp-icon-warn3[fill] {
  fill: var(--jp-warn-color3);
}

.jp-icon-warn0[stroke] {
  stroke: var(--jp-warn-color0);
}
.jp-icon-warn1[stroke] {
  stroke: var(--jp-warn-color1);
}
.jp-icon-warn2[stroke] {
  stroke: var(--jp-warn-color2);
}
.jp-icon-warn3[stroke] {
  stroke: var(--jp-warn-color3);
}
/* icon colors that contrast well with each other and most backgrounds */
.jp-icon-contrast0[fill] {
  fill: var(--jp-icon-contrast-color0);
}
.jp-icon-contrast1[fill] {
  fill: var(--jp-icon-contrast-color1);
}
.jp-icon-contrast2[fill] {
  fill: var(--jp-icon-contrast-color2);
}
.jp-icon-contrast3[fill] {
  fill: var(--jp-icon-contrast-color3);
}

.jp-icon-contrast0[stroke] {
  stroke: var(--jp-icon-contrast-color0);
}
.jp-icon-contrast1[stroke] {
  stroke: var(--jp-icon-contrast-color1);
}
.jp-icon-contrast2[stroke] {
  stroke: var(--jp-icon-contrast-color2);
}
.jp-icon-contrast3[stroke] {
  stroke: var(--jp-icon-contrast-color3);
}

/* CSS for icons in selected items in the settings editor */
#setting-editor .jp-PluginList .jp-mod-selected .jp-icon-selectable[fill] {
  fill: #fff;
}
#setting-editor
  .jp-PluginList
  .jp-mod-selected
  .jp-icon-selectable-inverse[fill] {
  fill: var(--jp-brand-color1);
}

/* CSS for icons in selected filebrowser listing items */
.jp-DirListing-item.jp-mod-selected .jp-icon-selectable[fill] {
  fill: #fff;
}
.jp-DirListing-item.jp-mod-selected .jp-icon-selectable-inverse[fill] {
  fill: var(--jp-brand-color1);
}

/* CSS for icons in selected tabs in the sidebar tab manager */
#tab-manager .lm-TabBar-tab.jp-mod-active .jp-icon-selectable[fill] {
  fill: #fff;
}

#tab-manager .lm-TabBar-tab.jp-mod-active .jp-icon-selectable-inverse[fill] {
  fill: var(--jp-brand-color1);
}
#tab-manager
  .lm-TabBar-tab.jp-mod-active
  .jp-icon-hover
  :hover
  .jp-icon-selectable[fill] {
  fill: var(--jp-brand-color1);
}

#tab-manager
  .lm-TabBar-tab.jp-mod-active
  .jp-icon-hover
  :hover
  .jp-icon-selectable-inverse[fill] {
  fill: #fff;
}

/**
 * TODO: come up with non css-hack solution for showing the busy icon on top
 *  of the close icon
 * CSS for complex behavior of close icon of tabs in the sidebar tab manager
 */
#tab-manager
  .lm-TabBar-tab.jp-mod-dirty
  > .lm-TabBar-tabCloseIcon
  > :not(:hover)
  > .jp-icon3[fill] {
  fill: none;
}
#tab-manager
  .lm-TabBar-tab.jp-mod-dirty
  > .lm-TabBar-tabCloseIcon
  > :not(:hover)
  > .jp-icon-busy[fill] {
  fill: var(--jp-inverse-layout-color3);
}

#tab-manager
  .lm-TabBar-tab.jp-mod-dirty.jp-mod-active
  > .lm-TabBar-tabCloseIcon
  > :not(:hover)
  > .jp-icon-busy[fill] {
  fill: #fff;
}

/**
* TODO: come up with non css-hack solution for showing the busy icon on top
*  of the close icon
* CSS for complex behavior of close icon of tabs in the main area tabbar
*/
.lm-DockPanel-tabBar
  .lm-TabBar-tab.lm-mod-closable.jp-mod-dirty
  > .lm-TabBar-tabCloseIcon
  > :not(:hover)
  > .jp-icon3[fill] {
  fill: none;
}
.lm-DockPanel-tabBar
  .lm-TabBar-tab.lm-mod-closable.jp-mod-dirty
  > .lm-TabBar-tabCloseIcon
  > :not(:hover)
  > .jp-icon-busy[fill] {
  fill: var(--jp-inverse-layout-color3);
}

/* CSS for icons in status bar */
#jp-main-statusbar .jp-mod-selected .jp-icon-selectable[fill] {
  fill: #fff;
}

#jp-main-statusbar .jp-mod-selected .jp-icon-selectable-inverse[fill] {
  fill: var(--jp-brand-color1);
}
/* special handling for splash icon CSS. While the theme CSS reloads during
   splash, the splash icon can loose theming. To prevent that, we set a
   default for its color variable */
:root {
  --jp-warn-color0: var(--md-orange-700);
}

/* not sure what to do with this one, used in filebrowser listing */
.jp-DragIcon {
  margin-right: 4px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/**
 * Support for alt colors for icons as inline SVG HTMLElements
 */

/* alt recolor the primary elements of an icon */
.jp-icon-alt .jp-icon0[fill] {
  fill: var(--jp-layout-color0);
}
.jp-icon-alt .jp-icon1[fill] {
  fill: var(--jp-layout-color1);
}
.jp-icon-alt .jp-icon2[fill] {
  fill: var(--jp-layout-color2);
}
.jp-icon-alt .jp-icon3[fill] {
  fill: var(--jp-layout-color3);
}
.jp-icon-alt .jp-icon4[fill] {
  fill: var(--jp-layout-color4);
}

.jp-icon-alt .jp-icon0[stroke] {
  stroke: var(--jp-layout-color0);
}
.jp-icon-alt .jp-icon1[stroke] {
  stroke: var(--jp-layout-color1);
}
.jp-icon-alt .jp-icon2[stroke] {
  stroke: var(--jp-layout-color2);
}
.jp-icon-alt .jp-icon3[stroke] {
  stroke: var(--jp-layout-color3);
}
.jp-icon-alt .jp-icon4[stroke] {
  stroke: var(--jp-layout-color4);
}

/* alt recolor the accent elements of an icon */
.jp-icon-alt .jp-icon-accent0[fill] {
  fill: var(--jp-inverse-layout-color0);
}
.jp-icon-alt .jp-icon-accent1[fill] {
  fill: var(--jp-inverse-layout-color1);
}
.jp-icon-alt .jp-icon-accent2[fill] {
  fill: var(--jp-inverse-layout-color2);
}
.jp-icon-alt .jp-icon-accent3[fill] {
  fill: var(--jp-inverse-layout-color3);
}
.jp-icon-alt .jp-icon-accent4[fill] {
  fill: var(--jp-inverse-layout-color4);
}

.jp-icon-alt .jp-icon-accent0[stroke] {
  stroke: var(--jp-inverse-layout-color0);
}
.jp-icon-alt .jp-icon-accent1[stroke] {
  stroke: var(--jp-inverse-layout-color1);
}
.jp-icon-alt .jp-icon-accent2[stroke] {
  stroke: var(--jp-inverse-layout-color2);
}
.jp-icon-alt .jp-icon-accent3[stroke] {
  stroke: var(--jp-inverse-layout-color3);
}
.jp-icon-alt .jp-icon-accent4[stroke] {
  stroke: var(--jp-inverse-layout-color4);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-icon-hoverShow:not(:hover) svg {
  display: none !important;
}

/**
 * Support for hover colors for icons as inline SVG HTMLElements
 */

/**
 * regular colors
 */

/* recolor the primary elements of an icon */
.jp-icon-hover :hover .jp-icon0-hover[fill] {
  fill: var(--jp-inverse-layout-color0);
}
.jp-icon-hover :hover .jp-icon1-hover[fill] {
  fill: var(--jp-inverse-layout-color1);
}
.jp-icon-hover :hover .jp-icon2-hover[fill] {
  fill: var(--jp-inverse-layout-color2);
}
.jp-icon-hover :hover .jp-icon3-hover[fill] {
  fill: var(--jp-inverse-layout-color3);
}
.jp-icon-hover :hover .jp-icon4-hover[fill] {
  fill: var(--jp-inverse-layout-color4);
}

.jp-icon-hover :hover .jp-icon0-hover[stroke] {
  stroke: var(--jp-inverse-layout-color0);
}
.jp-icon-hover :hover .jp-icon1-hover[stroke] {
  stroke: var(--jp-inverse-layout-color1);
}
.jp-icon-hover :hover .jp-icon2-hover[stroke] {
  stroke: var(--jp-inverse-layout-color2);
}
.jp-icon-hover :hover .jp-icon3-hover[stroke] {
  stroke: var(--jp-inverse-layout-color3);
}
.jp-icon-hover :hover .jp-icon4-hover[stroke] {
  stroke: var(--jp-inverse-layout-color4);
}

/* recolor the accent elements of an icon */
.jp-icon-hover :hover .jp-icon-accent0-hover[fill] {
  fill: var(--jp-layout-color0);
}
.jp-icon-hover :hover .jp-icon-accent1-hover[fill] {
  fill: var(--jp-layout-color1);
}
.jp-icon-hover :hover .jp-icon-accent2-hover[fill] {
  fill: var(--jp-layout-color2);
}
.jp-icon-hover :hover .jp-icon-accent3-hover[fill] {
  fill: var(--jp-layout-color3);
}
.jp-icon-hover :hover .jp-icon-accent4-hover[fill] {
  fill: var(--jp-layout-color4);
}

.jp-icon-hover :hover .jp-icon-accent0-hover[stroke] {
  stroke: var(--jp-layout-color0);
}
.jp-icon-hover :hover .jp-icon-accent1-hover[stroke] {
  stroke: var(--jp-layout-color1);
}
.jp-icon-hover :hover .jp-icon-accent2-hover[stroke] {
  stroke: var(--jp-layout-color2);
}
.jp-icon-hover :hover .jp-icon-accent3-hover[stroke] {
  stroke: var(--jp-layout-color3);
}
.jp-icon-hover :hover .jp-icon-accent4-hover[stroke] {
  stroke: var(--jp-layout-color4);
}

/* set the color of an icon to transparent */
.jp-icon-hover :hover .jp-icon-none-hover[fill] {
  fill: none;
}

.jp-icon-hover :hover .jp-icon-none-hover[stroke] {
  stroke: none;
}

/**
 * inverse colors
 */

/* inverse recolor the primary elements of an icon */
.jp-icon-hover.jp-icon-alt :hover .jp-icon0-hover[fill] {
  fill: var(--jp-layout-color0);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon1-hover[fill] {
  fill: var(--jp-layout-color1);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon2-hover[fill] {
  fill: var(--jp-layout-color2);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon3-hover[fill] {
  fill: var(--jp-layout-color3);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon4-hover[fill] {
  fill: var(--jp-layout-color4);
}

.jp-icon-hover.jp-icon-alt :hover .jp-icon0-hover[stroke] {
  stroke: var(--jp-layout-color0);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon1-hover[stroke] {
  stroke: var(--jp-layout-color1);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon2-hover[stroke] {
  stroke: var(--jp-layout-color2);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon3-hover[stroke] {
  stroke: var(--jp-layout-color3);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon4-hover[stroke] {
  stroke: var(--jp-layout-color4);
}

/* inverse recolor the accent elements of an icon */
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent0-hover[fill] {
  fill: var(--jp-inverse-layout-color0);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent1-hover[fill] {
  fill: var(--jp-inverse-layout-color1);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent2-hover[fill] {
  fill: var(--jp-inverse-layout-color2);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent3-hover[fill] {
  fill: var(--jp-inverse-layout-color3);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent4-hover[fill] {
  fill: var(--jp-inverse-layout-color4);
}

.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent0-hover[stroke] {
  stroke: var(--jp-inverse-layout-color0);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent1-hover[stroke] {
  stroke: var(--jp-inverse-layout-color1);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent2-hover[stroke] {
  stroke: var(--jp-inverse-layout-color2);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent3-hover[stroke] {
  stroke: var(--jp-inverse-layout-color3);
}
.jp-icon-hover.jp-icon-alt :hover .jp-icon-accent4-hover[stroke] {
  stroke: var(--jp-inverse-layout-color4);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* Sibling imports */

/* Override Blueprint's _reset.scss styles */
html {
  box-sizing: unset;
}

*,
*::before,
*::after {
  box-sizing: unset;
}

body {
  color: unset;
  font-family: var(--jp-ui-font-family);
}

p {
  margin-top: unset;
  margin-bottom: unset;
}

small {
  font-size: unset;
}

strong {
  font-weight: unset;
}

/* Override Blueprint's _typography.scss styles */
a {
  text-decoration: unset;
  color: unset;
}
a:hover {
  text-decoration: unset;
  color: unset;
}

/* Override Blueprint's _accessibility.scss styles */
:focus {
  outline: unset;
  outline-offset: unset;
  -moz-outline-radius: unset;
}

/* Styles for ui-components */
.jp-Button {
  border-radius: var(--jp-border-radius);
  padding: 0px 12px;
  font-size: var(--jp-ui-font-size1);
}

/* Use our own theme for hover styles */
button.jp-Button.bp3-button.bp3-minimal:hover {
  background-color: var(--jp-layout-color2);
}
.jp-Button.minimal {
  color: unset !important;
}

.jp-Button.jp-ToolbarButtonComponent {
  text-transform: none;
}

.jp-InputGroup input {
  box-sizing: border-box;
  border-radius: 0;
  background-color: transparent;
  color: var(--jp-ui-font-color0);
  box-shadow: inset 0 0 0 var(--jp-border-width) var(--jp-input-border-color);
}

.jp-InputGroup input:focus {
  box-shadow: inset 0 0 0 var(--jp-border-width)
      var(--jp-input-active-box-shadow-color),
    inset 0 0 0 3px var(--jp-input-active-box-shadow-color);
}

.jp-InputGroup input::placeholder,
input::placeholder {
  color: var(--jp-ui-font-color3);
}

.jp-BPIcon {
  display: inline-block;
  vertical-align: middle;
  margin: auto;
}

/* Stop blueprint futzing with our icon fills */
.bp3-icon.jp-BPIcon > svg:not([fill]) {
  fill: var(--jp-inverse-layout-color3);
}

.jp-InputGroupAction {
  padding: 6px;
}

.jp-HTMLSelect.jp-DefaultStyle select {
  background-color: initial;
  border: none;
  border-radius: 0;
  box-shadow: none;
  color: var(--jp-ui-font-color0);
  display: block;
  font-size: var(--jp-ui-font-size1);
  height: 24px;
  line-height: 14px;
  padding: 0 25px 0 10px;
  text-align: left;
  -moz-appearance: none;
  -webkit-appearance: none;
}

/* Use our own theme for hover and option styles */
.jp-HTMLSelect.jp-DefaultStyle select:hover,
.jp-HTMLSelect.jp-DefaultStyle select > option {
  background-color: var(--jp-layout-color2);
  color: var(--jp-ui-font-color0);
}
select {
  box-sizing: border-box;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-Collapse {
  display: flex;
  flex-direction: column;
  align-items: stretch;
  border-top: 1px solid var(--jp-border-color2);
  border-bottom: 1px solid var(--jp-border-color2);
}

.jp-Collapse-header {
  padding: 1px 12px;
  color: var(--jp-ui-font-color1);
  background-color: var(--jp-layout-color1);
  font-size: var(--jp-ui-font-size2);
}

.jp-Collapse-header:hover {
  background-color: var(--jp-layout-color2);
}

.jp-Collapse-contents {
  padding: 0px 12px 0px 12px;
  background-color: var(--jp-layout-color1);
  color: var(--jp-ui-font-color1);
  overflow: auto;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Variables
|----------------------------------------------------------------------------*/

:root {
  --jp-private-commandpalette-search-height: 28px;
}

/*-----------------------------------------------------------------------------
| Overall styles
|----------------------------------------------------------------------------*/

.lm-CommandPalette {
  padding-bottom: 0px;
  color: var(--jp-ui-font-color1);
  background: var(--jp-layout-color1);
  /* This is needed so that all font sizing of children done in ems is
   * relative to this base size */
  font-size: var(--jp-ui-font-size1);
}

/*-----------------------------------------------------------------------------
| Search
|----------------------------------------------------------------------------*/

.lm-CommandPalette-search {
  padding: 4px;
  background-color: var(--jp-layout-color1);
  z-index: 2;
}

.lm-CommandPalette-wrapper {
  overflow: overlay;
  padding: 0px 9px;
  background-color: var(--jp-input-active-background);
  height: 30px;
  box-shadow: inset 0 0 0 var(--jp-border-width) var(--jp-input-border-color);
}

.lm-CommandPalette.lm-mod-focused .lm-CommandPalette-wrapper {
  box-shadow: inset 0 0 0 1px var(--jp-input-active-box-shadow-color),
    inset 0 0 0 3px var(--jp-input-active-box-shadow-color);
}

.lm-CommandPalette-wrapper::after {
  content: ' ';
  color: white;
  background-color: var(--jp-brand-color1);
  position: absolute;
  top: 4px;
  right: 4px;
  height: 30px;
  width: 10px;
  padding: 0px 10px;
  background-image: var(--jp-icon-search-white);
  background-size: 20px;
  background-repeat: no-repeat;
  background-position: center;
}

.lm-CommandPalette-input {
  background: transparent;
  width: calc(100% - 18px);
  float: left;
  border: none;
  outline: none;
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color0);
  line-height: var(--jp-private-commandpalette-search-height);
}

.lm-CommandPalette-input::-webkit-input-placeholder,
.lm-CommandPalette-input::-moz-placeholder,
.lm-CommandPalette-input:-ms-input-placeholder {
  color: var(--jp-ui-font-color3);
  font-size: var(--jp-ui-font-size1);
}

/*-----------------------------------------------------------------------------
| Results
|----------------------------------------------------------------------------*/

.lm-CommandPalette-header:first-child {
  margin-top: 0px;
}

.lm-CommandPalette-header {
  border-bottom: solid var(--jp-border-width) var(--jp-border-color2);
  color: var(--jp-ui-font-color1);
  cursor: pointer;
  display: flex;
  font-size: var(--jp-ui-font-size0);
  font-weight: 600;
  letter-spacing: 1px;
  margin-top: 8px;
  padding: 8px 0 8px 12px;
  text-transform: uppercase;
}

.lm-CommandPalette-header.lm-mod-active {
  background: var(--jp-layout-color2);
}

.lm-CommandPalette-header > mark {
  background-color: transparent;
  font-weight: bold;
  color: var(--jp-ui-font-color1);
}

.lm-CommandPalette-item {
  padding: 4px 12px 4px 4px;
  color: var(--jp-ui-font-color1);
  font-size: var(--jp-ui-font-size1);
  font-weight: 400;
  display: flex;
}

.lm-CommandPalette-item.lm-mod-disabled {
  color: var(--jp-ui-font-color3);
}

.lm-CommandPalette-item.lm-mod-active {
  background: var(--jp-layout-color3);
}

.lm-CommandPalette-item.lm-mod-active:hover:not(.lm-mod-disabled) {
  background: var(--jp-layout-color4);
}

.lm-CommandPalette-item:hover:not(.lm-mod-active):not(.lm-mod-disabled) {
  background: var(--jp-layout-color2);
}

.lm-CommandPalette-itemContent {
  overflow: hidden;
}

.lm-CommandPalette-itemLabel > mark {
  color: var(--jp-ui-font-color0);
  background-color: transparent;
  font-weight: bold;
}

.lm-CommandPalette-item.lm-mod-disabled mark {
  color: var(--jp-ui-font-color3);
}

.lm-CommandPalette-item .lm-CommandPalette-itemIcon {
  margin: 0 4px 0 0;
  position: relative;
  width: 16px;
  top: 2px;
  flex: 0 0 auto;
}

.lm-CommandPalette-item.lm-mod-disabled .lm-CommandPalette-itemIcon {
  opacity: 0.4;
}

.lm-CommandPalette-item .lm-CommandPalette-itemShortcut {
  flex: 0 0 auto;
}

.lm-CommandPalette-itemCaption {
  display: none;
}

.lm-CommandPalette-content {
  background-color: var(--jp-layout-color1);
}

.lm-CommandPalette-content:empty:after {
  content: 'No results';
  margin: auto;
  margin-top: 20px;
  width: 100px;
  display: block;
  font-size: var(--jp-ui-font-size2);
  font-family: var(--jp-ui-font-family);
  font-weight: lighter;
}

.lm-CommandPalette-emptyMessage {
  text-align: center;
  margin-top: 24px;
  line-height: 1.32;
  padding: 0px 8px;
  color: var(--jp-content-font-color3);
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2017, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-Dialog {
  position: absolute;
  z-index: 10000;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  top: 0px;
  left: 0px;
  margin: 0;
  padding: 0;
  width: 100%;
  height: 100%;
  background: var(--jp-dialog-background);
}

.jp-Dialog-content {
  display: flex;
  flex-direction: column;
  margin-left: auto;
  margin-right: auto;
  background: var(--jp-layout-color1);
  padding: 24px;
  padding-bottom: 12px;
  min-width: 300px;
  min-height: 150px;
  max-width: 1000px;
  max-height: 500px;
  box-sizing: border-box;
  box-shadow: var(--jp-elevation-z20);
  word-wrap: break-word;
  border-radius: var(--jp-border-radius);
  /* This is needed so that all font sizing of children done in ems is
   * relative to this base size */
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color1);
}

.jp-Dialog-button {
  overflow: visible;
}

button.jp-Dialog-button:focus {
  outline: 1px solid var(--jp-brand-color1);
  outline-offset: 4px;
  -moz-outline-radius: 0px;
}

button.jp-Dialog-button:focus::-moz-focus-inner {
  border: 0;
}

.jp-Dialog-header {
  flex: 0 0 auto;
  padding-bottom: 12px;
  font-size: var(--jp-ui-font-size3);
  font-weight: 400;
  color: var(--jp-ui-font-color0);
}

.jp-Dialog-body {
  display: flex;
  flex-direction: column;
  flex: 1 1 auto;
  font-size: var(--jp-ui-font-size1);
  background: var(--jp-layout-color1);
  overflow: auto;
}

.jp-Dialog-footer {
  display: flex;
  flex-direction: row;
  justify-content: flex-end;
  flex: 0 0 auto;
  margin-left: -12px;
  margin-right: -12px;
  padding: 12px;
}

.jp-Dialog-title {
  overflow: hidden;
  white-space: nowrap;
  text-overflow: ellipsis;
}

.jp-Dialog-body > .jp-select-wrapper {
  width: 100%;
}

.jp-Dialog-body > button {
  padding: 0px 16px;
}

.jp-Dialog-body > label {
  line-height: 1.4;
  color: var(--jp-ui-font-color0);
}

.jp-Dialog-button.jp-mod-styled:not(:last-child) {
  margin-right: 12px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2016, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-HoverBox {
  position: fixed;
}

.jp-HoverBox.jp-mod-outofview {
  display: none;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-IFrame {
  width: 100%;
  height: 100%;
}

.jp-IFrame > iframe {
  border: none;
}

/*
When drag events occur, `p-mod-override-cursor` is added to the body.
Because iframes steal all cursor events, the following two rules are necessary
to suppress pointer events while resize drags are occurring. There may be a
better solution to this problem.
*/
body.lm-mod-override-cursor .jp-IFrame {
  position: relative;
}

body.lm-mod-override-cursor .jp-IFrame:before {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: transparent;
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2016, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-MainAreaWidget > :focus {
  outline: none;
}

/**
 * google-material-color v1.2.6
 * https://github.com/danlevan/google-material-color
 */
:root {
  --md-red-50: #ffebee;
  --md-red-100: #ffcdd2;
  --md-red-200: #ef9a9a;
  --md-red-300: #e57373;
  --md-red-400: #ef5350;
  --md-red-500: #f44336;
  --md-red-600: #e53935;
  --md-red-700: #d32f2f;
  --md-red-800: #c62828;
  --md-red-900: #b71c1c;
  --md-red-A100: #ff8a80;
  --md-red-A200: #ff5252;
  --md-red-A400: #ff1744;
  --md-red-A700: #d50000;

  --md-pink-50: #fce4ec;
  --md-pink-100: #f8bbd0;
  --md-pink-200: #f48fb1;
  --md-pink-300: #f06292;
  --md-pink-400: #ec407a;
  --md-pink-500: #e91e63;
  --md-pink-600: #d81b60;
  --md-pink-700: #c2185b;
  --md-pink-800: #ad1457;
  --md-pink-900: #880e4f;
  --md-pink-A100: #ff80ab;
  --md-pink-A200: #ff4081;
  --md-pink-A400: #f50057;
  --md-pink-A700: #c51162;

  --md-purple-50: #f3e5f5;
  --md-purple-100: #e1bee7;
  --md-purple-200: #ce93d8;
  --md-purple-300: #ba68c8;
  --md-purple-400: #ab47bc;
  --md-purple-500: #9c27b0;
  --md-purple-600: #8e24aa;
  --md-purple-700: #7b1fa2;
  --md-purple-800: #6a1b9a;
  --md-purple-900: #4a148c;
  --md-purple-A100: #ea80fc;
  --md-purple-A200: #e040fb;
  --md-purple-A400: #d500f9;
  --md-purple-A700: #aa00ff;

  --md-deep-purple-50: #ede7f6;
  --md-deep-purple-100: #d1c4e9;
  --md-deep-purple-200: #b39ddb;
  --md-deep-purple-300: #9575cd;
  --md-deep-purple-400: #7e57c2;
  --md-deep-purple-500: #673ab7;
  --md-deep-purple-600: #5e35b1;
  --md-deep-purple-700: #512da8;
  --md-deep-purple-800: #4527a0;
  --md-deep-purple-900: #311b92;
  --md-deep-purple-A100: #b388ff;
  --md-deep-purple-A200: #7c4dff;
  --md-deep-purple-A400: #651fff;
  --md-deep-purple-A700: #6200ea;

  --md-indigo-50: #e8eaf6;
  --md-indigo-100: #c5cae9;
  --md-indigo-200: #9fa8da;
  --md-indigo-300: #7986cb;
  --md-indigo-400: #5c6bc0;
  --md-indigo-500: #3f51b5;
  --md-indigo-600: #3949ab;
  --md-indigo-700: #303f9f;
  --md-indigo-800: #283593;
  --md-indigo-900: #1a237e;
  --md-indigo-A100: #8c9eff;
  --md-indigo-A200: #536dfe;
  --md-indigo-A400: #3d5afe;
  --md-indigo-A700: #304ffe;

  --md-blue-50: #e3f2fd;
  --md-blue-100: #bbdefb;
  --md-blue-200: #90caf9;
  --md-blue-300: #64b5f6;
  --md-blue-400: #42a5f5;
  --md-blue-500: #2196f3;
  --md-blue-600: #1e88e5;
  --md-blue-700: #1976d2;
  --md-blue-800: #1565c0;
  --md-blue-900: #0d47a1;
  --md-blue-A100: #82b1ff;
  --md-blue-A200: #448aff;
  --md-blue-A400: #2979ff;
  --md-blue-A700: #2962ff;

  --md-light-blue-50: #e1f5fe;
  --md-light-blue-100: #b3e5fc;
  --md-light-blue-200: #81d4fa;
  --md-light-blue-300: #4fc3f7;
  --md-light-blue-400: #29b6f6;
  --md-light-blue-500: #03a9f4;
  --md-light-blue-600: #039be5;
  --md-light-blue-700: #0288d1;
  --md-light-blue-800: #0277bd;
  --md-light-blue-900: #01579b;
  --md-light-blue-A100: #80d8ff;
  --md-light-blue-A200: #40c4ff;
  --md-light-blue-A400: #00b0ff;
  --md-light-blue-A700: #0091ea;

  --md-cyan-50: #e0f7fa;
  --md-cyan-100: #b2ebf2;
  --md-cyan-200: #80deea;
  --md-cyan-300: #4dd0e1;
  --md-cyan-400: #26c6da;
  --md-cyan-500: #00bcd4;
  --md-cyan-600: #00acc1;
  --md-cyan-700: #0097a7;
  --md-cyan-800: #00838f;
  --md-cyan-900: #006064;
  --md-cyan-A100: #84ffff;
  --md-cyan-A200: #18ffff;
  --md-cyan-A400: #00e5ff;
  --md-cyan-A700: #00b8d4;

  --md-teal-50: #e0f2f1;
  --md-teal-100: #b2dfdb;
  --md-teal-200: #80cbc4;
  --md-teal-300: #4db6ac;
  --md-teal-400: #26a69a;
  --md-teal-500: #009688;
  --md-teal-600: #00897b;
  --md-teal-700: #00796b;
  --md-teal-800: #00695c;
  --md-teal-900: #004d40;
  --md-teal-A100: #a7ffeb;
  --md-teal-A200: #64ffda;
  --md-teal-A400: #1de9b6;
  --md-teal-A700: #00bfa5;

  --md-green-50: #e8f5e9;
  --md-green-100: #c8e6c9;
  --md-green-200: #a5d6a7;
  --md-green-300: #81c784;
  --md-green-400: #66bb6a;
  --md-green-500: #4caf50;
  --md-green-600: #43a047;
  --md-green-700: #388e3c;
  --md-green-800: #2e7d32;
  --md-green-900: #1b5e20;
  --md-green-A100: #b9f6ca;
  --md-green-A200: #69f0ae;
  --md-green-A400: #00e676;
  --md-green-A700: #00c853;

  --md-light-green-50: #f1f8e9;
  --md-light-green-100: #dcedc8;
  --md-light-green-200: #c5e1a5;
  --md-light-green-300: #aed581;
  --md-light-green-400: #9ccc65;
  --md-light-green-500: #8bc34a;
  --md-light-green-600: #7cb342;
  --md-light-green-700: #689f38;
  --md-light-green-800: #558b2f;
  --md-light-green-900: #33691e;
  --md-light-green-A100: #ccff90;
  --md-light-green-A200: #b2ff59;
  --md-light-green-A400: #76ff03;
  --md-light-green-A700: #64dd17;

  --md-lime-50: #f9fbe7;
  --md-lime-100: #f0f4c3;
  --md-lime-200: #e6ee9c;
  --md-lime-300: #dce775;
  --md-lime-400: #d4e157;
  --md-lime-500: #cddc39;
  --md-lime-600: #c0ca33;
  --md-lime-700: #afb42b;
  --md-lime-800: #9e9d24;
  --md-lime-900: #827717;
  --md-lime-A100: #f4ff81;
  --md-lime-A200: #eeff41;
  --md-lime-A400: #c6ff00;
  --md-lime-A700: #aeea00;

  --md-yellow-50: #fffde7;
  --md-yellow-100: #fff9c4;
  --md-yellow-200: #fff59d;
  --md-yellow-300: #fff176;
  --md-yellow-400: #ffee58;
  --md-yellow-500: #ffeb3b;
  --md-yellow-600: #fdd835;
  --md-yellow-700: #fbc02d;
  --md-yellow-800: #f9a825;
  --md-yellow-900: #f57f17;
  --md-yellow-A100: #ffff8d;
  --md-yellow-A200: #ffff00;
  --md-yellow-A400: #ffea00;
  --md-yellow-A700: #ffd600;

  --md-amber-50: #fff8e1;
  --md-amber-100: #ffecb3;
  --md-amber-200: #ffe082;
  --md-amber-300: #ffd54f;
  --md-amber-400: #ffca28;
  --md-amber-500: #ffc107;
  --md-amber-600: #ffb300;
  --md-amber-700: #ffa000;
  --md-amber-800: #ff8f00;
  --md-amber-900: #ff6f00;
  --md-amber-A100: #ffe57f;
  --md-amber-A200: #ffd740;
  --md-amber-A400: #ffc400;
  --md-amber-A700: #ffab00;

  --md-orange-50: #fff3e0;
  --md-orange-100: #ffe0b2;
  --md-orange-200: #ffcc80;
  --md-orange-300: #ffb74d;
  --md-orange-400: #ffa726;
  --md-orange-500: #ff9800;
  --md-orange-600: #fb8c00;
  --md-orange-700: #f57c00;
  --md-orange-800: #ef6c00;
  --md-orange-900: #e65100;
  --md-orange-A100: #ffd180;
  --md-orange-A200: #ffab40;
  --md-orange-A400: #ff9100;
  --md-orange-A700: #ff6d00;

  --md-deep-orange-50: #fbe9e7;
  --md-deep-orange-100: #ffccbc;
  --md-deep-orange-200: #ffab91;
  --md-deep-orange-300: #ff8a65;
  --md-deep-orange-400: #ff7043;
  --md-deep-orange-500: #ff5722;
  --md-deep-orange-600: #f4511e;
  --md-deep-orange-700: #e64a19;
  --md-deep-orange-800: #d84315;
  --md-deep-orange-900: #bf360c;
  --md-deep-orange-A100: #ff9e80;
  --md-deep-orange-A200: #ff6e40;
  --md-deep-orange-A400: #ff3d00;
  --md-deep-orange-A700: #dd2c00;

  --md-brown-50: #efebe9;
  --md-brown-100: #d7ccc8;
  --md-brown-200: #bcaaa4;
  --md-brown-300: #a1887f;
  --md-brown-400: #8d6e63;
  --md-brown-500: #795548;
  --md-brown-600: #6d4c41;
  --md-brown-700: #5d4037;
  --md-brown-800: #4e342e;
  --md-brown-900: #3e2723;

  --md-grey-50: #fafafa;
  --md-grey-100: #f5f5f5;
  --md-grey-200: #eeeeee;
  --md-grey-300: #e0e0e0;
  --md-grey-400: #bdbdbd;
  --md-grey-500: #9e9e9e;
  --md-grey-600: #757575;
  --md-grey-700: #616161;
  --md-grey-800: #424242;
  --md-grey-900: #212121;

  --md-blue-grey-50: #eceff1;
  --md-blue-grey-100: #cfd8dc;
  --md-blue-grey-200: #b0bec5;
  --md-blue-grey-300: #90a4ae;
  --md-blue-grey-400: #78909c;
  --md-blue-grey-500: #607d8b;
  --md-blue-grey-600: #546e7a;
  --md-blue-grey-700: #455a64;
  --md-blue-grey-800: #37474f;
  --md-blue-grey-900: #263238;
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2017, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-Spinner {
  position: absolute;
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 10;
  left: 0;
  top: 0;
  width: 100%;
  height: 100%;
  background: var(--jp-layout-color0);
  outline: none;
}

.jp-SpinnerContent {
  font-size: 10px;
  margin: 50px auto;
  text-indent: -9999em;
  width: 3em;
  height: 3em;
  border-radius: 50%;
  background: var(--jp-brand-color3);
  background: linear-gradient(
    to right,
    #f37626 10%,
    rgba(255, 255, 255, 0) 42%
  );
  position: relative;
  animation: load3 1s infinite linear, fadeIn 1s;
}

.jp-SpinnerContent:before {
  width: 50%;
  height: 50%;
  background: #f37626;
  border-radius: 100% 0 0 0;
  position: absolute;
  top: 0;
  left: 0;
  content: '';
}

.jp-SpinnerContent:after {
  background: var(--jp-layout-color0);
  width: 75%;
  height: 75%;
  border-radius: 50%;
  content: '';
  margin: auto;
  position: absolute;
  top: 0;
  left: 0;
  bottom: 0;
  right: 0;
}

@keyframes fadeIn {
  0% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

@keyframes load3 {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
  }
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2017, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

button.jp-mod-styled {
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color0);
  border: none;
  box-sizing: border-box;
  text-align: center;
  line-height: 32px;
  height: 32px;
  padding: 0px 12px;
  letter-spacing: 0.8px;
  outline: none;
  appearance: none;
  -webkit-appearance: none;
  -moz-appearance: none;
}

input.jp-mod-styled {
  background: var(--jp-input-background);
  height: 28px;
  box-sizing: border-box;
  border: var(--jp-border-width) solid var(--jp-border-color1);
  padding-left: 7px;
  padding-right: 7px;
  font-size: var(--jp-ui-font-size2);
  color: var(--jp-ui-font-color0);
  outline: none;
  appearance: none;
  -webkit-appearance: none;
  -moz-appearance: none;
}

input.jp-mod-styled:focus {
  border: var(--jp-border-width) solid var(--md-blue-500);
  box-shadow: inset 0 0 4px var(--md-blue-300);
}

.jp-select-wrapper {
  display: flex;
  position: relative;
  flex-direction: column;
  padding: 1px;
  background-color: var(--jp-layout-color1);
  height: 28px;
  box-sizing: border-box;
  margin-bottom: 12px;
}

.jp-select-wrapper.jp-mod-focused select.jp-mod-styled {
  border: var(--jp-border-width) solid var(--jp-input-active-border-color);
  box-shadow: var(--jp-input-box-shadow);
  background-color: var(--jp-input-active-background);
}

select.jp-mod-styled:hover {
  background-color: var(--jp-layout-color1);
  cursor: pointer;
  color: var(--jp-ui-font-color0);
  background-color: var(--jp-input-hover-background);
  box-shadow: inset 0 0px 1px rgba(0, 0, 0, 0.5);
}

select.jp-mod-styled {
  flex: 1 1 auto;
  height: 32px;
  width: 100%;
  font-size: var(--jp-ui-font-size2);
  background: var(--jp-input-background);
  color: var(--jp-ui-font-color0);
  padding: 0 25px 0 8px;
  border: var(--jp-border-width) solid var(--jp-input-border-color);
  border-radius: 0px;
  outline: none;
  appearance: none;
  -webkit-appearance: none;
  -moz-appearance: none;
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2016, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

:root {
  --jp-private-toolbar-height: calc(
    28px + var(--jp-border-width)
  ); /* leave 28px for content */
}

.jp-Toolbar {
  color: var(--jp-ui-font-color1);
  flex: 0 0 auto;
  display: flex;
  flex-direction: row;
  border-bottom: var(--jp-border-width) solid var(--jp-toolbar-border-color);
  box-shadow: var(--jp-toolbar-box-shadow);
  background: var(--jp-toolbar-background);
  min-height: var(--jp-toolbar-micro-height);
  padding: 2px;
  z-index: 1;
}

/* Toolbar items */

.jp-Toolbar > .jp-Toolbar-item.jp-Toolbar-spacer {
  flex-grow: 1;
  flex-shrink: 1;
}

.jp-Toolbar-item.jp-Toolbar-kernelStatus {
  display: inline-block;
  width: 32px;
  background-repeat: no-repeat;
  background-position: center;
  background-size: 16px;
}

.jp-Toolbar > .jp-Toolbar-item {
  flex: 0 0 auto;
  display: flex;
  padding-left: 1px;
  padding-right: 1px;
  font-size: var(--jp-ui-font-size1);
  line-height: var(--jp-private-toolbar-height);
  height: 100%;
}

/* Toolbar buttons */

/* This is the div we use to wrap the react component into a Widget */
div.jp-ToolbarButton {
  color: transparent;
  border: none;
  box-sizing: border-box;
  outline: none;
  appearance: none;
  -webkit-appearance: none;
  -moz-appearance: none;
  padding: 0px;
  margin: 0px;
}

button.jp-ToolbarButtonComponent {
  background: var(--jp-layout-color1);
  border: none;
  box-sizing: border-box;
  outline: none;
  appearance: none;
  -webkit-appearance: none;
  -moz-appearance: none;
  padding: 0px 6px;
  margin: 0px;
  height: 24px;
  border-radius: var(--jp-border-radius);
  display: flex;
  align-items: center;
  text-align: center;
  font-size: 14px;
  min-width: unset;
  min-height: unset;
}

button.jp-ToolbarButtonComponent:disabled {
  opacity: 0.4;
}

button.jp-ToolbarButtonComponent span {
  padding: 0px;
  flex: 0 0 auto;
}

button.jp-ToolbarButtonComponent .jp-ToolbarButtonComponent-label {
  font-size: var(--jp-ui-font-size1);
  line-height: 100%;
  padding-left: 2px;
  color: var(--jp-ui-font-color1);
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2017, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Copyright (c) 2014-2017, PhosphorJS Contributors
|
| Distributed under the terms of the BSD 3-Clause License.
|
| The full license is in the file LICENSE, distributed with this software.
|----------------------------------------------------------------------------*/


/* <DEPRECATED> */ body.p-mod-override-cursor *, /* </DEPRECATED> */
body.lm-mod-override-cursor * {
  cursor: inherit !important;
}

/*-----------------------------------------------------------------------------
| Copyright (c) 2014-2016, Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-JSONEditor {
  display: flex;
  flex-direction: column;
  width: 100%;
}

.jp-JSONEditor-host {
  flex: 1 1 auto;
  border: var(--jp-border-width) solid var(--jp-input-border-color);
  border-radius: 0px;
  background: var(--jp-layout-color0);
  min-height: 50px;
  padding: 1px;
}

.jp-JSONEditor.jp-mod-error .jp-JSONEditor-host {
  border-color: red;
  outline-color: red;
}

.jp-JSONEditor-header {
  display: flex;
  flex: 1 0 auto;
  padding: 0 0 0 12px;
}

.jp-JSONEditor-header label {
  flex: 0 0 auto;
}

.jp-JSONEditor-commitButton {
  height: 16px;
  width: 16px;
  background-size: 18px;
  background-repeat: no-repeat;
  background-position: center;
}

.jp-JSONEditor-host.jp-mod-focused {
  background-color: var(--jp-input-active-background);
  border: 1px solid var(--jp-input-active-border-color);
  box-shadow: var(--jp-input-box-shadow);
}

.jp-Editor.jp-mod-dropTarget {
  border: var(--jp-border-width) solid var(--jp-input-active-border-color);
  box-shadow: var(--jp-input-box-shadow);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* BASICS */

.CodeMirror {
  /* Set height, width, borders, and global font properties here */
  font-family: monospace;
  height: 300px;
  color: black;
  direction: ltr;
}

/* PADDING */

.CodeMirror-lines {
  padding: 4px 0; /* Vertical padding around content */
}
.CodeMirror pre.CodeMirror-line,
.CodeMirror pre.CodeMirror-line-like {
  padding: 0 4px; /* Horizontal padding of content */
}

.CodeMirror-scrollbar-filler, .CodeMirror-gutter-filler {
  background-color: white; /* The little square between H and V scrollbars */
}

/* GUTTER */

.CodeMirror-gutters {
  border-right: 1px solid #ddd;
  background-color: #f7f7f7;
  white-space: nowrap;
}
.CodeMirror-linenumbers {}
.CodeMirror-linenumber {
  padding: 0 3px 0 5px;
  min-width: 20px;
  text-align: right;
  color: #999;
  white-space: nowrap;
}

.CodeMirror-guttermarker { color: black; }
.CodeMirror-guttermarker-subtle { color: #999; }

/* CURSOR */

.CodeMirror-cursor {
  border-left: 1px solid black;
  border-right: none;
  width: 0;
}
/* Shown when moving in bi-directional text */
.CodeMirror div.CodeMirror-secondarycursor {
  border-left: 1px solid silver;
}
.cm-fat-cursor .CodeMirror-cursor {
  width: auto;
  border: 0 !important;
  background: #7e7;
}
.cm-fat-cursor div.CodeMirror-cursors {
  z-index: 1;
}
.cm-fat-cursor-mark {
  background-color: rgba(20, 255, 20, 0.5);
  -webkit-animation: blink 1.06s steps(1) infinite;
  -moz-animation: blink 1.06s steps(1) infinite;
  animation: blink 1.06s steps(1) infinite;
}
.cm-animate-fat-cursor {
  width: auto;
  border: 0;
  -webkit-animation: blink 1.06s steps(1) infinite;
  -moz-animation: blink 1.06s steps(1) infinite;
  animation: blink 1.06s steps(1) infinite;
  background-color: #7e7;
}
@-moz-keyframes blink {
  0% {}
  50% { background-color: transparent; }
  100% {}
}
@-webkit-keyframes blink {
  0% {}
  50% { background-color: transparent; }
  100% {}
}
@keyframes blink {
  0% {}
  50% { background-color: transparent; }
  100% {}
}

/* Can style cursor different in overwrite (non-insert) mode */
.CodeMirror-overwrite .CodeMirror-cursor {}

.cm-tab { display: inline-block; text-decoration: inherit; }

.CodeMirror-rulers {
  position: absolute;
  left: 0; right: 0; top: -50px; bottom: 0;
  overflow: hidden;
}
.CodeMirror-ruler {
  border-left: 1px solid #ccc;
  top: 0; bottom: 0;
  position: absolute;
}

/* DEFAULT THEME */

.cm-s-default .cm-header {color: blue;}
.cm-s-default .cm-quote {color: #090;}
.cm-negative {color: #d44;}
.cm-positive {color: #292;}
.cm-header, .cm-strong {font-weight: bold;}
.cm-em {font-style: italic;}
.cm-link {text-decoration: underline;}
.cm-strikethrough {text-decoration: line-through;}

.cm-s-default .cm-keyword {color: #708;}
.cm-s-default .cm-atom {color: #219;}
.cm-s-default .cm-number {color: #164;}
.cm-s-default .cm-def {color: #00f;}
.cm-s-default .cm-variable,
.cm-s-default .cm-punctuation,
.cm-s-default .cm-property,
.cm-s-default .cm-operator {}
.cm-s-default .cm-variable-2 {color: #05a;}
.cm-s-default .cm-variable-3, .cm-s-default .cm-type {color: #085;}
.cm-s-default .cm-comment {color: #a50;}
.cm-s-default .cm-string {color: #a11;}
.cm-s-default .cm-string-2 {color: #f50;}
.cm-s-default .cm-meta {color: #555;}
.cm-s-default .cm-qualifier {color: #555;}
.cm-s-default .cm-builtin {color: #30a;}
.cm-s-default .cm-bracket {color: #997;}
.cm-s-default .cm-tag {color: #170;}
.cm-s-default .cm-attribute {color: #00c;}
.cm-s-default .cm-hr {color: #999;}
.cm-s-default .cm-link {color: #00c;}

.cm-s-default .cm-error {color: #f00;}
.cm-invalidchar {color: #f00;}

.CodeMirror-composing { border-bottom: 2px solid; }

/* Default styles for common addons */

div.CodeMirror span.CodeMirror-matchingbracket {color: #0b0;}
div.CodeMirror span.CodeMirror-nonmatchingbracket {color: #a22;}
.CodeMirror-matchingtag { background: rgba(255, 150, 0, .3); }
.CodeMirror-activeline-background {background: #e8f2ff;}

/* STOP */

/* The rest of this file contains styles related to the mechanics of
   the editor. You probably shouldn't touch them. */

.CodeMirror {
  position: relative;
  overflow: hidden;
  background: white;
}

.CodeMirror-scroll {
  overflow: scroll !important; /* Things will break if this is overridden */
  /* 30px is the magic margin used to hide the element's real scrollbars */
  /* See overflow: hidden in .CodeMirror */
  margin-bottom: -30px; margin-right: -30px;
  padding-bottom: 30px;
  height: 100%;
  outline: none; /* Prevent dragging from highlighting the element */
  position: relative;
}
.CodeMirror-sizer {
  position: relative;
  border-right: 30px solid transparent;
}

/* The fake, visible scrollbars. Used to force redraw during scrolling
   before actual scrolling happens, thus preventing shaking and
   flickering artifacts. */
.CodeMirror-vscrollbar, .CodeMirror-hscrollbar, .CodeMirror-scrollbar-filler, .CodeMirror-gutter-filler {
  position: absolute;
  z-index: 6;
  display: none;
}
.CodeMirror-vscrollbar {
  right: 0; top: 0;
  overflow-x: hidden;
  overflow-y: scroll;
}
.CodeMirror-hscrollbar {
  bottom: 0; left: 0;
  overflow-y: hidden;
  overflow-x: scroll;
}
.CodeMirror-scrollbar-filler {
  right: 0; bottom: 0;
}
.CodeMirror-gutter-filler {
  left: 0; bottom: 0;
}

.CodeMirror-gutters {
  position: absolute; left: 0; top: 0;
  min-height: 100%;
  z-index: 3;
}
.CodeMirror-gutter {
  white-space: normal;
  height: 100%;
  display: inline-block;
  vertical-align: top;
  margin-bottom: -30px;
}
.CodeMirror-gutter-wrapper {
  position: absolute;
  z-index: 4;
  background: none !important;
  border: none !important;
}
.CodeMirror-gutter-background {
  position: absolute;
  top: 0; bottom: 0;
  z-index: 4;
}
.CodeMirror-gutter-elt {
  position: absolute;
  cursor: default;
  z-index: 4;
}
.CodeMirror-gutter-wrapper ::selection { background-color: transparent }
.CodeMirror-gutter-wrapper ::-moz-selection { background-color: transparent }

.CodeMirror-lines {
  cursor: text;
  min-height: 1px; /* prevents collapsing before first draw */
}
.CodeMirror pre.CodeMirror-line,
.CodeMirror pre.CodeMirror-line-like {
  /* Reset some styles that the rest of the page might have set */
  -moz-border-radius: 0; -webkit-border-radius: 0; border-radius: 0;
  border-width: 0;
  background: transparent;
  font-family: inherit;
  font-size: inherit;
  margin: 0;
  white-space: pre;
  word-wrap: normal;
  line-height: inherit;
  color: inherit;
  z-index: 2;
  position: relative;
  overflow: visible;
  -webkit-tap-highlight-color: transparent;
  -webkit-font-variant-ligatures: contextual;
  font-variant-ligatures: contextual;
}
.CodeMirror-wrap pre.CodeMirror-line,
.CodeMirror-wrap pre.CodeMirror-line-like {
  word-wrap: break-word;
  white-space: pre-wrap;
  word-break: normal;
}

.CodeMirror-linebackground {
  position: absolute;
  left: 0; right: 0; top: 0; bottom: 0;
  z-index: 0;
}

.CodeMirror-linewidget {
  position: relative;
  z-index: 2;
  padding: 0.1px; /* Force widget margins to stay inside of the container */
}

.CodeMirror-widget {}

.CodeMirror-rtl pre { direction: rtl; }

.CodeMirror-code {
  outline: none;
}

/* Force content-box sizing for the elements where we expect it */
.CodeMirror-scroll,
.CodeMirror-sizer,
.CodeMirror-gutter,
.CodeMirror-gutters,
.CodeMirror-linenumber {
  -moz-box-sizing: content-box;
  box-sizing: content-box;
}

.CodeMirror-measure {
  position: absolute;
  width: 100%;
  height: 0;
  overflow: hidden;
  visibility: hidden;
}

.CodeMirror-cursor {
  position: absolute;
  pointer-events: none;
}
.CodeMirror-measure pre { position: static; }

div.CodeMirror-cursors {
  visibility: hidden;
  position: relative;
  z-index: 3;
}
div.CodeMirror-dragcursors {
  visibility: visible;
}

.CodeMirror-focused div.CodeMirror-cursors {
  visibility: visible;
}

.CodeMirror-selected { background: #d9d9d9; }
.CodeMirror-focused .CodeMirror-selected { background: #d7d4f0; }
.CodeMirror-crosshair { cursor: crosshair; }
.CodeMirror-line::selection, .CodeMirror-line > span::selection, .CodeMirror-line > span > span::selection { background: #d7d4f0; }
.CodeMirror-line::-moz-selection, .CodeMirror-line > span::-moz-selection, .CodeMirror-line > span > span::-moz-selection { background: #d7d4f0; }

.cm-searching {
  background-color: #ffa;
  background-color: rgba(255, 255, 0, .4);
}

/* Used to force a border model for a node */
.cm-force-border { padding-right: .1px; }

@media print {
  /* Hide the cursor when printing */
  .CodeMirror div.CodeMirror-cursors {
    visibility: hidden;
  }
}

/* See issue #2901 */
.cm-tab-wrap-hack:after { content: ''; }

/* Help users use markselection to safely style text background */
span.CodeMirror-selectedtext { background: none; }

.CodeMirror-dialog {
  position: absolute;
  left: 0; right: 0;
  background: inherit;
  z-index: 15;
  padding: .1em .8em;
  overflow: hidden;
  color: inherit;
}

.CodeMirror-dialog-top {
  border-bottom: 1px solid #eee;
  top: 0;
}

.CodeMirror-dialog-bottom {
  border-top: 1px solid #eee;
  bottom: 0;
}

.CodeMirror-dialog input {
  border: none;
  outline: none;
  background: transparent;
  width: 20em;
  color: inherit;
  font-family: monospace;
}

.CodeMirror-dialog button {
  font-size: 70%;
}

.CodeMirror-foldmarker {
  color: blue;
  text-shadow: #b9f 1px 1px 2px, #b9f -1px -1px 2px, #b9f 1px -1px 2px, #b9f -1px 1px 2px;
  font-family: arial;
  line-height: .3;
  cursor: pointer;
}
.CodeMirror-foldgutter {
  width: .7em;
}
.CodeMirror-foldgutter-open,
.CodeMirror-foldgutter-folded {
  cursor: pointer;
}
.CodeMirror-foldgutter-open:after {
  content: "\25BE";
}
.CodeMirror-foldgutter-folded:after {
  content: "\25B8";
}

/*
  Name:       material
  Author:     Mattia Astorino (http://github.com/equinusocio)
  Website:    https://material-theme.site/
*/

.cm-s-material.CodeMirror {
  background-color: #263238;
  color: #EEFFFF;
}

.cm-s-material .CodeMirror-gutters {
  background: #263238;
  color: #546E7A;
  border: none;
}

.cm-s-material .CodeMirror-guttermarker,
.cm-s-material .CodeMirror-guttermarker-subtle,
.cm-s-material .CodeMirror-linenumber {
  color: #546E7A;
}

.cm-s-material .CodeMirror-cursor {
  border-left: 1px solid #FFCC00;
}

.cm-s-material div.CodeMirror-selected {
  background: rgba(128, 203, 196, 0.2);
}

.cm-s-material.CodeMirror-focused div.CodeMirror-selected {
  background: rgba(128, 203, 196, 0.2);
}

.cm-s-material .CodeMirror-line::selection,
.cm-s-material .CodeMirror-line>span::selection,
.cm-s-material .CodeMirror-line>span>span::selection {
  background: rgba(128, 203, 196, 0.2);
}

.cm-s-material .CodeMirror-line::-moz-selection,
.cm-s-material .CodeMirror-line>span::-moz-selection,
.cm-s-material .CodeMirror-line>span>span::-moz-selection {
  background: rgba(128, 203, 196, 0.2);
}

.cm-s-material .CodeMirror-activeline-background {
  background: rgba(0, 0, 0, 0.5);
}

.cm-s-material .cm-keyword {
  color: #C792EA;
}

.cm-s-material .cm-operator {
  color: #89DDFF;
}

.cm-s-material .cm-variable-2 {
  color: #EEFFFF;
}

.cm-s-material .cm-variable-3,
.cm-s-material .cm-type {
  color: #f07178;
}

.cm-s-material .cm-builtin {
  color: #FFCB6B;
}

.cm-s-material .cm-atom {
  color: #F78C6C;
}

.cm-s-material .cm-number {
  color: #FF5370;
}

.cm-s-material .cm-def {
  color: #82AAFF;
}

.cm-s-material .cm-string {
  color: #C3E88D;
}

.cm-s-material .cm-string-2 {
  color: #f07178;
}

.cm-s-material .cm-comment {
  color: #546E7A;
}

.cm-s-material .cm-variable {
  color: #f07178;
}

.cm-s-material .cm-tag {
  color: #FF5370;
}

.cm-s-material .cm-meta {
  color: #FFCB6B;
}

.cm-s-material .cm-attribute {
  color: #C792EA;
}

.cm-s-material .cm-property {
  color: #C792EA;
}

.cm-s-material .cm-qualifier {
  color: #DECB6B;
}

.cm-s-material .cm-variable-3,
.cm-s-material .cm-type {
  color: #DECB6B;
}


.cm-s-material .cm-error {
  color: rgba(255, 255, 255, 1.0);
  background-color: #FF5370;
}

.cm-s-material .CodeMirror-matchingbracket {
  text-decoration: underline;
  color: white !important;
}
/**
 * "
 *  Using Zenburn color palette from the Emacs Zenburn Theme
 *  https://github.com/bbatsov/zenburn-emacs/blob/master/zenburn-theme.el
 *
 *  Also using parts of https://github.com/xavi/coderay-lighttable-theme
 * "
 * From: https://github.com/wisenomad/zenburn-lighttable-theme/blob/master/zenburn.css
 */

.cm-s-zenburn .CodeMirror-gutters { background: #3f3f3f !important; }
.cm-s-zenburn .CodeMirror-foldgutter-open, .CodeMirror-foldgutter-folded { color: #999; }
.cm-s-zenburn .CodeMirror-cursor { border-left: 1px solid white; }
.cm-s-zenburn { background-color: #3f3f3f; color: #dcdccc; }
.cm-s-zenburn span.cm-builtin { color: #dcdccc; font-weight: bold; }
.cm-s-zenburn span.cm-comment { color: #7f9f7f; }
.cm-s-zenburn span.cm-keyword { color: #f0dfaf; font-weight: bold; }
.cm-s-zenburn span.cm-atom { color: #bfebbf; }
.cm-s-zenburn span.cm-def { color: #dcdccc; }
.cm-s-zenburn span.cm-variable { color: #dfaf8f; }
.cm-s-zenburn span.cm-variable-2 { color: #dcdccc; }
.cm-s-zenburn span.cm-string { color: #cc9393; }
.cm-s-zenburn span.cm-string-2 { color: #cc9393; }
.cm-s-zenburn span.cm-number { color: #dcdccc; }
.cm-s-zenburn span.cm-tag { color: #93e0e3; }
.cm-s-zenburn span.cm-property { color: #dfaf8f; }
.cm-s-zenburn span.cm-attribute { color: #dfaf8f; }
.cm-s-zenburn span.cm-qualifier { color: #7cb8bb; }
.cm-s-zenburn span.cm-meta { color: #f0dfaf; }
.cm-s-zenburn span.cm-header { color: #f0efd0; }
.cm-s-zenburn span.cm-operator { color: #f0efd0; }
.cm-s-zenburn span.CodeMirror-matchingbracket { box-sizing: border-box; background: transparent; border-bottom: 1px solid; }
.cm-s-zenburn span.CodeMirror-nonmatchingbracket { border-bottom: 1px solid; background: none; }
.cm-s-zenburn .CodeMirror-activeline { background: #000000; }
.cm-s-zenburn .CodeMirror-activeline-background { background: #000000; }
.cm-s-zenburn div.CodeMirror-selected { background: #545454; }
.cm-s-zenburn .CodeMirror-focused div.CodeMirror-selected { background: #4f4f4f; }

.cm-s-abcdef.CodeMirror { background: #0f0f0f; color: #defdef; }
.cm-s-abcdef div.CodeMirror-selected { background: #515151; }
.cm-s-abcdef .CodeMirror-line::selection, .cm-s-abcdef .CodeMirror-line > span::selection, .cm-s-abcdef .CodeMirror-line > span > span::selection { background: rgba(56, 56, 56, 0.99); }
.cm-s-abcdef .CodeMirror-line::-moz-selection, .cm-s-abcdef .CodeMirror-line > span::-moz-selection, .cm-s-abcdef .CodeMirror-line > span > span::-moz-selection { background: rgba(56, 56, 56, 0.99); }
.cm-s-abcdef .CodeMirror-gutters { background: #555; border-right: 2px solid #314151; }
.cm-s-abcdef .CodeMirror-guttermarker { color: #222; }
.cm-s-abcdef .CodeMirror-guttermarker-subtle { color: azure; }
.cm-s-abcdef .CodeMirror-linenumber { color: #FFFFFF; }
.cm-s-abcdef .CodeMirror-cursor { border-left: 1px solid #00FF00; }

.cm-s-abcdef span.cm-keyword { color: darkgoldenrod; font-weight: bold; }
.cm-s-abcdef span.cm-atom { color: #77F; }
.cm-s-abcdef span.cm-number { color: violet; }
.cm-s-abcdef span.cm-def { color: #fffabc; }
.cm-s-abcdef span.cm-variable { color: #abcdef; }
.cm-s-abcdef span.cm-variable-2 { color: #cacbcc; }
.cm-s-abcdef span.cm-variable-3, .cm-s-abcdef span.cm-type { color: #def; }
.cm-s-abcdef span.cm-property { color: #fedcba; }
.cm-s-abcdef span.cm-operator { color: #ff0; }
.cm-s-abcdef span.cm-comment { color: #7a7b7c; font-style: italic;}
.cm-s-abcdef span.cm-string { color: #2b4; }
.cm-s-abcdef span.cm-meta { color: #C9F; }
.cm-s-abcdef span.cm-qualifier { color: #FFF700; }
.cm-s-abcdef span.cm-builtin { color: #30aabc; }
.cm-s-abcdef span.cm-bracket { color: #8a8a8a; }
.cm-s-abcdef span.cm-tag { color: #FFDD44; }
.cm-s-abcdef span.cm-attribute { color: #DDFF00; }
.cm-s-abcdef span.cm-error { color: #FF0000; }
.cm-s-abcdef span.cm-header { color: aquamarine; font-weight: bold; }
.cm-s-abcdef span.cm-link { color: blueviolet; }

.cm-s-abcdef .CodeMirror-activeline-background { background: #314151; }

/*

    Name:       Base16 Default Light
    Author:     Chris Kempson (http://chriskempson.com)

    CodeMirror template by Jan T. Sott (https://github.com/idleberg/base16-codemirror)
    Original Base16 color scheme by Chris Kempson (https://github.com/chriskempson/base16)

*/

.cm-s-base16-light.CodeMirror { background: #f5f5f5; color: #202020; }
.cm-s-base16-light div.CodeMirror-selected { background: #e0e0e0; }
.cm-s-base16-light .CodeMirror-line::selection, .cm-s-base16-light .CodeMirror-line > span::selection, .cm-s-base16-light .CodeMirror-line > span > span::selection { background: #e0e0e0; }
.cm-s-base16-light .CodeMirror-line::-moz-selection, .cm-s-base16-light .CodeMirror-line > span::-moz-selection, .cm-s-base16-light .CodeMirror-line > span > span::-moz-selection { background: #e0e0e0; }
.cm-s-base16-light .CodeMirror-gutters { background: #f5f5f5; border-right: 0px; }
.cm-s-base16-light .CodeMirror-guttermarker { color: #ac4142; }
.cm-s-base16-light .CodeMirror-guttermarker-subtle { color: #b0b0b0; }
.cm-s-base16-light .CodeMirror-linenumber { color: #b0b0b0; }
.cm-s-base16-light .CodeMirror-cursor { border-left: 1px solid #505050; }

.cm-s-base16-light span.cm-comment { color: #8f5536; }
.cm-s-base16-light span.cm-atom { color: #aa759f; }
.cm-s-base16-light span.cm-number { color: #aa759f; }

.cm-s-base16-light span.cm-property, .cm-s-base16-light span.cm-attribute { color: #90a959; }
.cm-s-base16-light span.cm-keyword { color: #ac4142; }
.cm-s-base16-light span.cm-string { color: #f4bf75; }

.cm-s-base16-light span.cm-variable { color: #90a959; }
.cm-s-base16-light span.cm-variable-2 { color: #6a9fb5; }
.cm-s-base16-light span.cm-def { color: #d28445; }
.cm-s-base16-light span.cm-bracket { color: #202020; }
.cm-s-base16-light span.cm-tag { color: #ac4142; }
.cm-s-base16-light span.cm-link { color: #aa759f; }
.cm-s-base16-light span.cm-error { background: #ac4142; color: #505050; }

.cm-s-base16-light .CodeMirror-activeline-background { background: #DDDCDC; }
.cm-s-base16-light .CodeMirror-matchingbracket { color: #f5f5f5 !important; background-color: #6A9FB5 !important}

/*

    Name:       Base16 Default Dark
    Author:     Chris Kempson (http://chriskempson.com)

    CodeMirror template by Jan T. Sott (https://github.com/idleberg/base16-codemirror)
    Original Base16 color scheme by Chris Kempson (https://github.com/chriskempson/base16)

*/

.cm-s-base16-dark.CodeMirror { background: #151515; color: #e0e0e0; }
.cm-s-base16-dark div.CodeMirror-selected { background: #303030; }
.cm-s-base16-dark .CodeMirror-line::selection, .cm-s-base16-dark .CodeMirror-line > span::selection, .cm-s-base16-dark .CodeMirror-line > span > span::selection { background: rgba(48, 48, 48, .99); }
.cm-s-base16-dark .CodeMirror-line::-moz-selection, .cm-s-base16-dark .CodeMirror-line > span::-moz-selection, .cm-s-base16-dark .CodeMirror-line > span > span::-moz-selection { background: rgba(48, 48, 48, .99); }
.cm-s-base16-dark .CodeMirror-gutters { background: #151515; border-right: 0px; }
.cm-s-base16-dark .CodeMirror-guttermarker { color: #ac4142; }
.cm-s-base16-dark .CodeMirror-guttermarker-subtle { color: #505050; }
.cm-s-base16-dark .CodeMirror-linenumber { color: #505050; }
.cm-s-base16-dark .CodeMirror-cursor { border-left: 1px solid #b0b0b0; }

.cm-s-base16-dark span.cm-comment { color: #8f5536; }
.cm-s-base16-dark span.cm-atom { color: #aa759f; }
.cm-s-base16-dark span.cm-number { color: #aa759f; }

.cm-s-base16-dark span.cm-property, .cm-s-base16-dark span.cm-attribute { color: #90a959; }
.cm-s-base16-dark span.cm-keyword { color: #ac4142; }
.cm-s-base16-dark span.cm-string { color: #f4bf75; }

.cm-s-base16-dark span.cm-variable { color: #90a959; }
.cm-s-base16-dark span.cm-variable-2 { color: #6a9fb5; }
.cm-s-base16-dark span.cm-def { color: #d28445; }
.cm-s-base16-dark span.cm-bracket { color: #e0e0e0; }
.cm-s-base16-dark span.cm-tag { color: #ac4142; }
.cm-s-base16-dark span.cm-link { color: #aa759f; }
.cm-s-base16-dark span.cm-error { background: #ac4142; color: #b0b0b0; }

.cm-s-base16-dark .CodeMirror-activeline-background { background: #202020; }
.cm-s-base16-dark .CodeMirror-matchingbracket { text-decoration: underline; color: white !important; }

/*

    Name:       dracula
    Author:     Michael Kaminsky (http://github.com/mkaminsky11)

    Original dracula color scheme by Zeno Rocha (https://github.com/zenorocha/dracula-theme)

*/


.cm-s-dracula.CodeMirror, .cm-s-dracula .CodeMirror-gutters {
  background-color: #282a36 !important;
  color: #f8f8f2 !important;
  border: none;
}
.cm-s-dracula .CodeMirror-gutters { color: #282a36; }
.cm-s-dracula .CodeMirror-cursor { border-left: solid thin #f8f8f0; }
.cm-s-dracula .CodeMirror-linenumber { color: #6D8A88; }
.cm-s-dracula .CodeMirror-selected { background: rgba(255, 255, 255, 0.10); }
.cm-s-dracula .CodeMirror-line::selection, .cm-s-dracula .CodeMirror-line > span::selection, .cm-s-dracula .CodeMirror-line > span > span::selection { background: rgba(255, 255, 255, 0.10); }
.cm-s-dracula .CodeMirror-line::-moz-selection, .cm-s-dracula .CodeMirror-line > span::-moz-selection, .cm-s-dracula .CodeMirror-line > span > span::-moz-selection { background: rgba(255, 255, 255, 0.10); }
.cm-s-dracula span.cm-comment { color: #6272a4; }
.cm-s-dracula span.cm-string, .cm-s-dracula span.cm-string-2 { color: #f1fa8c; }
.cm-s-dracula span.cm-number { color: #bd93f9; }
.cm-s-dracula span.cm-variable { color: #50fa7b; }
.cm-s-dracula span.cm-variable-2 { color: white; }
.cm-s-dracula span.cm-def { color: #50fa7b; }
.cm-s-dracula span.cm-operator { color: #ff79c6; }
.cm-s-dracula span.cm-keyword { color: #ff79c6; }
.cm-s-dracula span.cm-atom { color: #bd93f9; }
.cm-s-dracula span.cm-meta { color: #f8f8f2; }
.cm-s-dracula span.cm-tag { color: #ff79c6; }
.cm-s-dracula span.cm-attribute { color: #50fa7b; }
.cm-s-dracula span.cm-qualifier { color: #50fa7b; }
.cm-s-dracula span.cm-property { color: #66d9ef; }
.cm-s-dracula span.cm-builtin { color: #50fa7b; }
.cm-s-dracula span.cm-variable-3, .cm-s-dracula span.cm-type { color: #ffb86c; }

.cm-s-dracula .CodeMirror-activeline-background { background: rgba(255,255,255,0.1); }
.cm-s-dracula .CodeMirror-matchingbracket { text-decoration: underline; color: white !important; }

/*

    Name:       Hopscotch
    Author:     Jan T. Sott

    CodeMirror template by Jan T. Sott (https://github.com/idleberg/base16-codemirror)
    Original Base16 color scheme by Chris Kempson (https://github.com/chriskempson/base16)

*/

.cm-s-hopscotch.CodeMirror {background: #322931; color: #d5d3d5;}
.cm-s-hopscotch div.CodeMirror-selected {background: #433b42 !important;}
.cm-s-hopscotch .CodeMirror-gutters {background: #322931; border-right: 0px;}
.cm-s-hopscotch .CodeMirror-linenumber {color: #797379;}
.cm-s-hopscotch .CodeMirror-cursor {border-left: 1px solid #989498 !important;}

.cm-s-hopscotch span.cm-comment {color: #b33508;}
.cm-s-hopscotch span.cm-atom {color: #c85e7c;}
.cm-s-hopscotch span.cm-number {color: #c85e7c;}

.cm-s-hopscotch span.cm-property, .cm-s-hopscotch span.cm-attribute {color: #8fc13e;}
.cm-s-hopscotch span.cm-keyword {color: #dd464c;}
.cm-s-hopscotch span.cm-string {color: #fdcc59;}

.cm-s-hopscotch span.cm-variable {color: #8fc13e;}
.cm-s-hopscotch span.cm-variable-2 {color: #1290bf;}
.cm-s-hopscotch span.cm-def {color: #fd8b19;}
.cm-s-hopscotch span.cm-error {background: #dd464c; color: #989498;}
.cm-s-hopscotch span.cm-bracket {color: #d5d3d5;}
.cm-s-hopscotch span.cm-tag {color: #dd464c;}
.cm-s-hopscotch span.cm-link {color: #c85e7c;}

.cm-s-hopscotch .CodeMirror-matchingbracket { text-decoration: underline; color: white !important;}
.cm-s-hopscotch .CodeMirror-activeline-background { background: #302020; }

/****************************************************************/
/*   Based on mbonaci's Brackets mbo theme                      */
/*   https://github.com/mbonaci/global/blob/master/Mbo.tmTheme  */
/*   Create your own: http://tmtheme-editor.herokuapp.com       */
/****************************************************************/

.cm-s-mbo.CodeMirror { background: #2c2c2c; color: #ffffec; }
.cm-s-mbo div.CodeMirror-selected { background: #716C62; }
.cm-s-mbo .CodeMirror-line::selection, .cm-s-mbo .CodeMirror-line > span::selection, .cm-s-mbo .CodeMirror-line > span > span::selection { background: rgba(113, 108, 98, .99); }
.cm-s-mbo .CodeMirror-line::-moz-selection, .cm-s-mbo .CodeMirror-line > span::-moz-selection, .cm-s-mbo .CodeMirror-line > span > span::-moz-selection { background: rgba(113, 108, 98, .99); }
.cm-s-mbo .CodeMirror-gutters { background: #4e4e4e; border-right: 0px; }
.cm-s-mbo .CodeMirror-guttermarker { color: white; }
.cm-s-mbo .CodeMirror-guttermarker-subtle { color: grey; }
.cm-s-mbo .CodeMirror-linenumber { color: #dadada; }
.cm-s-mbo .CodeMirror-cursor { border-left: 1px solid #ffffec; }

.cm-s-mbo span.cm-comment { color: #95958a; }
.cm-s-mbo span.cm-atom { color: #00a8c6; }
.cm-s-mbo span.cm-number { color: #00a8c6; }

.cm-s-mbo span.cm-property, .cm-s-mbo span.cm-attribute { color: #9ddfe9; }
.cm-s-mbo span.cm-keyword { color: #ffb928; }
.cm-s-mbo span.cm-string { color: #ffcf6c; }
.cm-s-mbo span.cm-string.cm-property { color: #ffffec; }

.cm-s-mbo span.cm-variable { color: #ffffec; }
.cm-s-mbo span.cm-variable-2 { color: #00a8c6; }
.cm-s-mbo span.cm-def { color: #ffffec; }
.cm-s-mbo span.cm-bracket { color: #fffffc; font-weight: bold; }
.cm-s-mbo span.cm-tag { color: #9ddfe9; }
.cm-s-mbo span.cm-link { color: #f54b07; }
.cm-s-mbo span.cm-error { border-bottom: #636363; color: #ffffec; }
.cm-s-mbo span.cm-qualifier { color: #ffffec; }

.cm-s-mbo .CodeMirror-activeline-background { background: #494b41; }
.cm-s-mbo .CodeMirror-matchingbracket { color: #ffb928 !important; }
.cm-s-mbo .CodeMirror-matchingtag { background: rgba(255, 255, 255, .37); }

/*
  MDN-LIKE Theme - Mozilla
  Ported to CodeMirror by Peter Kroon <plakroon@gmail.com>
  Report bugs/issues here: https://github.com/codemirror/CodeMirror/issues
  GitHub: @peterkroon

  The mdn-like theme is inspired on the displayed code examples at: https://developer.mozilla.org/en-US/docs/Web/CSS/animation

*/
.cm-s-mdn-like.CodeMirror { color: #999; background-color: #fff; }
.cm-s-mdn-like div.CodeMirror-selected { background: #cfc; }
.cm-s-mdn-like .CodeMirror-line::selection, .cm-s-mdn-like .CodeMirror-line > span::selection, .cm-s-mdn-like .CodeMirror-line > span > span::selection { background: #cfc; }
.cm-s-mdn-like .CodeMirror-line::-moz-selection, .cm-s-mdn-like .CodeMirror-line > span::-moz-selection, .cm-s-mdn-like .CodeMirror-line > span > span::-moz-selection { background: #cfc; }

.cm-s-mdn-like .CodeMirror-gutters { background: #f8f8f8; border-left: 6px solid rgba(0,83,159,0.65); color: #333; }
.cm-s-mdn-like .CodeMirror-linenumber { color: #aaa; padding-left: 8px; }
.cm-s-mdn-like .CodeMirror-cursor { border-left: 2px solid #222; }

.cm-s-mdn-like .cm-keyword { color: #6262FF; }
.cm-s-mdn-like .cm-atom { color: #F90; }
.cm-s-mdn-like .cm-number { color:  #ca7841; }
.cm-s-mdn-like .cm-def { color: #8DA6CE; }
.cm-s-mdn-like span.cm-variable-2, .cm-s-mdn-like span.cm-tag { color: #690; }
.cm-s-mdn-like span.cm-variable-3, .cm-s-mdn-like span.cm-def, .cm-s-mdn-like span.cm-type { color: #07a; }

.cm-s-mdn-like .cm-variable { color: #07a; }
.cm-s-mdn-like .cm-property { color: #905; }
.cm-s-mdn-like .cm-qualifier { color: #690; }

.cm-s-mdn-like .cm-operator { color: #cda869; }
.cm-s-mdn-like .cm-comment { color:#777; font-weight:normal; }
.cm-s-mdn-like .cm-string { color:#07a; font-style:italic; }
.cm-s-mdn-like .cm-string-2 { color:#bd6b18; } /*?*/
.cm-s-mdn-like .cm-meta { color: #000; } /*?*/
.cm-s-mdn-like .cm-builtin { color: #9B7536; } /*?*/
.cm-s-mdn-like .cm-tag { color: #997643; }
.cm-s-mdn-like .cm-attribute { color: #d6bb6d; } /*?*/
.cm-s-mdn-like .cm-header { color: #FF6400; }
.cm-s-mdn-like .cm-hr { color: #AEAEAE; }
.cm-s-mdn-like .cm-link { color:#ad9361; font-style:italic; text-decoration:none; }
.cm-s-mdn-like .cm-error { border-bottom: 1px solid red; }

div.cm-s-mdn-like .CodeMirror-activeline-background { background: #efefff; }
div.cm-s-mdn-like span.CodeMirror-matchingbracket { outline:1px solid grey; color: inherit; }

.cm-s-mdn-like.CodeMirror { background-image: url(data:image/png;base64,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); }

/*

    Name:       seti
    Author:     Michael Kaminsky (http://github.com/mkaminsky11)

    Original seti color scheme by Jesse Weed (https://github.com/jesseweed/seti-syntax)

*/


.cm-s-seti.CodeMirror {
  background-color: #151718 !important;
  color: #CFD2D1 !important;
  border: none;
}
.cm-s-seti .CodeMirror-gutters {
  color: #404b53;
  background-color: #0E1112;
  border: none;
}
.cm-s-seti .CodeMirror-cursor { border-left: solid thin #f8f8f0; }
.cm-s-seti .CodeMirror-linenumber { color: #6D8A88; }
.cm-s-seti.CodeMirror-focused div.CodeMirror-selected { background: rgba(255, 255, 255, 0.10); }
.cm-s-seti .CodeMirror-line::selection, .cm-s-seti .CodeMirror-line > span::selection, .cm-s-seti .CodeMirror-line > span > span::selection { background: rgba(255, 255, 255, 0.10); }
.cm-s-seti .CodeMirror-line::-moz-selection, .cm-s-seti .CodeMirror-line > span::-moz-selection, .cm-s-seti .CodeMirror-line > span > span::-moz-selection { background: rgba(255, 255, 255, 0.10); }
.cm-s-seti span.cm-comment { color: #41535b; }
.cm-s-seti span.cm-string, .cm-s-seti span.cm-string-2 { color: #55b5db; }
.cm-s-seti span.cm-number { color: #cd3f45; }
.cm-s-seti span.cm-variable { color: #55b5db; }
.cm-s-seti span.cm-variable-2 { color: #a074c4; }
.cm-s-seti span.cm-def { color: #55b5db; }
.cm-s-seti span.cm-keyword { color: #ff79c6; }
.cm-s-seti span.cm-operator { color: #9fca56; }
.cm-s-seti span.cm-keyword { color: #e6cd69; }
.cm-s-seti span.cm-atom { color: #cd3f45; }
.cm-s-seti span.cm-meta { color: #55b5db; }
.cm-s-seti span.cm-tag { color: #55b5db; }
.cm-s-seti span.cm-attribute { color: #9fca56; }
.cm-s-seti span.cm-qualifier { color: #9fca56; }
.cm-s-seti span.cm-property { color: #a074c4; }
.cm-s-seti span.cm-variable-3, .cm-s-seti span.cm-type { color: #9fca56; }
.cm-s-seti span.cm-builtin { color: #9fca56; }
.cm-s-seti .CodeMirror-activeline-background { background: #101213; }
.cm-s-seti .CodeMirror-matchingbracket { text-decoration: underline; color: white !important; }

/*
Solarized theme for code-mirror
http://ethanschoonover.com/solarized
*/

/*
Solarized color palette
http://ethanschoonover.com/solarized/img/solarized-palette.png
*/

.solarized.base03 { color: #002b36; }
.solarized.base02 { color: #073642; }
.solarized.base01 { color: #586e75; }
.solarized.base00 { color: #657b83; }
.solarized.base0 { color: #839496; }
.solarized.base1 { color: #93a1a1; }
.solarized.base2 { color: #eee8d5; }
.solarized.base3  { color: #fdf6e3; }
.solarized.solar-yellow  { color: #b58900; }
.solarized.solar-orange  { color: #cb4b16; }
.solarized.solar-red { color: #dc322f; }
.solarized.solar-magenta { color: #d33682; }
.solarized.solar-violet  { color: #6c71c4; }
.solarized.solar-blue { color: #268bd2; }
.solarized.solar-cyan { color: #2aa198; }
.solarized.solar-green { color: #859900; }

/* Color scheme for code-mirror */

.cm-s-solarized {
  line-height: 1.45em;
  color-profile: sRGB;
  rendering-intent: auto;
}
.cm-s-solarized.cm-s-dark {
  color: #839496;
  background-color: #002b36;
  text-shadow: #002b36 0 1px;
}
.cm-s-solarized.cm-s-light {
  background-color: #fdf6e3;
  color: #657b83;
  text-shadow: #eee8d5 0 1px;
}

.cm-s-solarized .CodeMirror-widget {
  text-shadow: none;
}

.cm-s-solarized .cm-header { color: #586e75; }
.cm-s-solarized .cm-quote { color: #93a1a1; }

.cm-s-solarized .cm-keyword { color: #cb4b16; }
.cm-s-solarized .cm-atom { color: #d33682; }
.cm-s-solarized .cm-number { color: #d33682; }
.cm-s-solarized .cm-def { color: #2aa198; }

.cm-s-solarized .cm-variable { color: #839496; }
.cm-s-solarized .cm-variable-2 { color: #b58900; }
.cm-s-solarized .cm-variable-3, .cm-s-solarized .cm-type { color: #6c71c4; }

.cm-s-solarized .cm-property { color: #2aa198; }
.cm-s-solarized .cm-operator { color: #6c71c4; }

.cm-s-solarized .cm-comment { color: #586e75; font-style:italic; }

.cm-s-solarized .cm-string { color: #859900; }
.cm-s-solarized .cm-string-2 { color: #b58900; }

.cm-s-solarized .cm-meta { color: #859900; }
.cm-s-solarized .cm-qualifier { color: #b58900; }
.cm-s-solarized .cm-builtin { color: #d33682; }
.cm-s-solarized .cm-bracket { color: #cb4b16; }
.cm-s-solarized .CodeMirror-matchingbracket { color: #859900; }
.cm-s-solarized .CodeMirror-nonmatchingbracket { color: #dc322f; }
.cm-s-solarized .cm-tag { color: #93a1a1; }
.cm-s-solarized .cm-attribute { color: #2aa198; }
.cm-s-solarized .cm-hr {
  color: transparent;
  border-top: 1px solid #586e75;
  display: block;
}
.cm-s-solarized .cm-link { color: #93a1a1; cursor: pointer; }
.cm-s-solarized .cm-special { color: #6c71c4; }
.cm-s-solarized .cm-em {
  color: #999;
  text-decoration: underline;
  text-decoration-style: dotted;
}
.cm-s-solarized .cm-error,
.cm-s-solarized .cm-invalidchar {
  color: #586e75;
  border-bottom: 1px dotted #dc322f;
}

.cm-s-solarized.cm-s-dark div.CodeMirror-selected { background: #073642; }
.cm-s-solarized.cm-s-dark.CodeMirror ::selection { background: rgba(7, 54, 66, 0.99); }
.cm-s-solarized.cm-s-dark .CodeMirror-line::-moz-selection, .cm-s-dark .CodeMirror-line > span::-moz-selection, .cm-s-dark .CodeMirror-line > span > span::-moz-selection { background: rgba(7, 54, 66, 0.99); }

.cm-s-solarized.cm-s-light div.CodeMirror-selected { background: #eee8d5; }
.cm-s-solarized.cm-s-light .CodeMirror-line::selection, .cm-s-light .CodeMirror-line > span::selection, .cm-s-light .CodeMirror-line > span > span::selection { background: #eee8d5; }
.cm-s-solarized.cm-s-light .CodeMirror-line::-moz-selection, .cm-s-ligh .CodeMirror-line > span::-moz-selection, .cm-s-ligh .CodeMirror-line > span > span::-moz-selection { background: #eee8d5; }

/* Editor styling */



/* Little shadow on the view-port of the buffer view */
.cm-s-solarized.CodeMirror {
  -moz-box-shadow: inset 7px 0 12px -6px #000;
  -webkit-box-shadow: inset 7px 0 12px -6px #000;
  box-shadow: inset 7px 0 12px -6px #000;
}

/* Remove gutter border */
.cm-s-solarized .CodeMirror-gutters {
  border-right: 0;
}

/* Gutter colors and line number styling based of color scheme (dark / light) */

/* Dark */
.cm-s-solarized.cm-s-dark .CodeMirror-gutters {
  background-color: #073642;
}

.cm-s-solarized.cm-s-dark .CodeMirror-linenumber {
  color: #586e75;
  text-shadow: #021014 0 -1px;
}

/* Light */
.cm-s-solarized.cm-s-light .CodeMirror-gutters {
  background-color: #eee8d5;
}

.cm-s-solarized.cm-s-light .CodeMirror-linenumber {
  color: #839496;
}

/* Common */
.cm-s-solarized .CodeMirror-linenumber {
  padding: 0 5px;
}
.cm-s-solarized .CodeMirror-guttermarker-subtle { color: #586e75; }
.cm-s-solarized.cm-s-dark .CodeMirror-guttermarker { color: #ddd; }
.cm-s-solarized.cm-s-light .CodeMirror-guttermarker { color: #cb4b16; }

.cm-s-solarized .CodeMirror-gutter .CodeMirror-gutter-text {
  color: #586e75;
}

/* Cursor */
.cm-s-solarized .CodeMirror-cursor { border-left: 1px solid #819090; }

/* Fat cursor */
.cm-s-solarized.cm-s-light.cm-fat-cursor .CodeMirror-cursor { background: #77ee77; }
.cm-s-solarized.cm-s-light .cm-animate-fat-cursor { background-color: #77ee77; }
.cm-s-solarized.cm-s-dark.cm-fat-cursor .CodeMirror-cursor { background: #586e75; }
.cm-s-solarized.cm-s-dark .cm-animate-fat-cursor { background-color: #586e75; }

/* Active line */
.cm-s-solarized.cm-s-dark .CodeMirror-activeline-background {
  background: rgba(255, 255, 255, 0.06);
}
.cm-s-solarized.cm-s-light .CodeMirror-activeline-background {
  background: rgba(0, 0, 0, 0.06);
}

.cm-s-the-matrix.CodeMirror { background: #000000; color: #00FF00; }
.cm-s-the-matrix div.CodeMirror-selected { background: #2D2D2D; }
.cm-s-the-matrix .CodeMirror-line::selection, .cm-s-the-matrix .CodeMirror-line > span::selection, .cm-s-the-matrix .CodeMirror-line > span > span::selection { background: rgba(45, 45, 45, 0.99); }
.cm-s-the-matrix .CodeMirror-line::-moz-selection, .cm-s-the-matrix .CodeMirror-line > span::-moz-selection, .cm-s-the-matrix .CodeMirror-line > span > span::-moz-selection { background: rgba(45, 45, 45, 0.99); }
.cm-s-the-matrix .CodeMirror-gutters { background: #060; border-right: 2px solid #00FF00; }
.cm-s-the-matrix .CodeMirror-guttermarker { color: #0f0; }
.cm-s-the-matrix .CodeMirror-guttermarker-subtle { color: white; }
.cm-s-the-matrix .CodeMirror-linenumber { color: #FFFFFF; }
.cm-s-the-matrix .CodeMirror-cursor { border-left: 1px solid #00FF00; }

.cm-s-the-matrix span.cm-keyword { color: #008803; font-weight: bold; }
.cm-s-the-matrix span.cm-atom { color: #3FF; }
.cm-s-the-matrix span.cm-number { color: #FFB94F; }
.cm-s-the-matrix span.cm-def { color: #99C; }
.cm-s-the-matrix span.cm-variable { color: #F6C; }
.cm-s-the-matrix span.cm-variable-2 { color: #C6F; }
.cm-s-the-matrix span.cm-variable-3, .cm-s-the-matrix span.cm-type { color: #96F; }
.cm-s-the-matrix span.cm-property { color: #62FFA0; }
.cm-s-the-matrix span.cm-operator { color: #999; }
.cm-s-the-matrix span.cm-comment { color: #CCCCCC; }
.cm-s-the-matrix span.cm-string { color: #39C; }
.cm-s-the-matrix span.cm-meta { color: #C9F; }
.cm-s-the-matrix span.cm-qualifier { color: #FFF700; }
.cm-s-the-matrix span.cm-builtin { color: #30a; }
.cm-s-the-matrix span.cm-bracket { color: #cc7; }
.cm-s-the-matrix span.cm-tag { color: #FFBD40; }
.cm-s-the-matrix span.cm-attribute { color: #FFF700; }
.cm-s-the-matrix span.cm-error { color: #FF0000; }

.cm-s-the-matrix .CodeMirror-activeline-background { background: #040; }

/*
Copyright (C) 2011 by MarkLogic Corporation
Author: Mike Brevoort <mike@brevoort.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
.cm-s-xq-light span.cm-keyword { line-height: 1em; font-weight: bold; color: #5A5CAD; }
.cm-s-xq-light span.cm-atom { color: #6C8CD5; }
.cm-s-xq-light span.cm-number { color: #164; }
.cm-s-xq-light span.cm-def { text-decoration:underline; }
.cm-s-xq-light span.cm-variable { color: black; }
.cm-s-xq-light span.cm-variable-2 { color:black; }
.cm-s-xq-light span.cm-variable-3, .cm-s-xq-light span.cm-type { color: black; }
.cm-s-xq-light span.cm-property {}
.cm-s-xq-light span.cm-operator {}
.cm-s-xq-light span.cm-comment { color: #0080FF; font-style: italic; }
.cm-s-xq-light span.cm-string { color: red; }
.cm-s-xq-light span.cm-meta { color: yellow; }
.cm-s-xq-light span.cm-qualifier { color: grey; }
.cm-s-xq-light span.cm-builtin { color: #7EA656; }
.cm-s-xq-light span.cm-bracket { color: #cc7; }
.cm-s-xq-light span.cm-tag { color: #3F7F7F; }
.cm-s-xq-light span.cm-attribute { color: #7F007F; }
.cm-s-xq-light span.cm-error { color: #f00; }

.cm-s-xq-light .CodeMirror-activeline-background { background: #e8f2ff; }
.cm-s-xq-light .CodeMirror-matchingbracket { outline:1px solid grey;color:black !important;background:yellow; }

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.CodeMirror {
  line-height: var(--jp-code-line-height);
  font-size: var(--jp-code-font-size);
  font-family: var(--jp-code-font-family);
  border: 0;
  border-radius: 0;
  height: auto;
  /* Changed to auto to autogrow */
}

.CodeMirror pre {
  padding: 0 var(--jp-code-padding);
}

.jp-CodeMirrorEditor[data-type='inline'] .CodeMirror-dialog {
  background-color: var(--jp-layout-color0);
  color: var(--jp-content-font-color1);
}

/* This causes https://github.com/jupyter/jupyterlab/issues/522 */
/* May not cause it not because we changed it! */
.CodeMirror-lines {
  padding: var(--jp-code-padding) 0;
}

.CodeMirror-linenumber {
  padding: 0 8px;
}

.jp-CodeMirrorEditor-static {
  margin: var(--jp-code-padding);
}

.jp-CodeMirrorEditor,
.jp-CodeMirrorEditor-static {
  cursor: text;
}

.jp-CodeMirrorEditor[data-type='inline'] .CodeMirror-cursor {
  border-left: var(--jp-code-cursor-width0) solid var(--jp-editor-cursor-color);
}

/* When zoomed out 67% and 33% on a screen of 1440 width x 900 height */
@media screen and (min-width: 2138px) and (max-width: 4319px) {
  .jp-CodeMirrorEditor[data-type='inline'] .CodeMirror-cursor {
    border-left: var(--jp-code-cursor-width1) solid
      var(--jp-editor-cursor-color);
  }
}

/* When zoomed out less than 33% */
@media screen and (min-width: 4320px) {
  .jp-CodeMirrorEditor[data-type='inline'] .CodeMirror-cursor {
    border-left: var(--jp-code-cursor-width2) solid
      var(--jp-editor-cursor-color);
  }
}

.CodeMirror.jp-mod-readOnly .CodeMirror-cursor {
  display: none;
}

.CodeMirror-gutters {
  border-right: 1px solid var(--jp-border-color2);
  background-color: var(--jp-layout-color0);
}

.jp-CollaboratorCursor {
  border-left: 5px solid transparent;
  border-right: 5px solid transparent;
  border-top: none;
  border-bottom: 3px solid;
  background-clip: content-box;
  margin-left: -5px;
  margin-right: -5px;
}

.CodeMirror-selectedtext.cm-searching {
  background-color: var(--jp-search-selected-match-background-color) !important;
  color: var(--jp-search-selected-match-color) !important;
}

.cm-searching {
  background-color: var(
    --jp-search-unselected-match-background-color
  ) !important;
  color: var(--jp-search-unselected-match-color) !important;
}

.CodeMirror-focused .CodeMirror-selected {
  background-color: var(--jp-editor-selected-focused-background);
}

.CodeMirror-selected {
  background-color: var(--jp-editor-selected-background);
}

.jp-CollaboratorCursor-hover {
  position: absolute;
  z-index: 1;
  transform: translateX(-50%);
  color: white;
  border-radius: 3px;
  padding-left: 4px;
  padding-right: 4px;
  padding-top: 1px;
  padding-bottom: 1px;
  text-align: center;
  font-size: var(--jp-ui-font-size1);
  white-space: nowrap;
}

.jp-CodeMirror-ruler {
  border-left: 1px dashed var(--jp-border-color2);
}

/**
 * Here is our jupyter theme for CodeMirror syntax highlighting
 * This is used in our marked.js syntax highlighting and CodeMirror itself
 * The string "jupyter" is set in ../codemirror/widget.DEFAULT_CODEMIRROR_THEME
 * This came from the classic notebook, which came form highlight.js/GitHub
 */

/**
 * CodeMirror themes are handling the background/color in this way. This works
 * fine for CodeMirror editors outside the notebook, but the notebook styles
 * these things differently.
 */
.CodeMirror.cm-s-jupyter {
  background: var(--jp-layout-color0);
  color: var(--jp-content-font-color1);
}

/* In the notebook, we want this styling to be handled by its container */
.jp-CodeConsole .CodeMirror.cm-s-jupyter,
.jp-Notebook .CodeMirror.cm-s-jupyter {
  background: transparent;
}

.cm-s-jupyter .CodeMirror-cursor {
  border-left: var(--jp-code-cursor-width0) solid var(--jp-editor-cursor-color);
}
.cm-s-jupyter span.cm-keyword {
  color: var(--jp-mirror-editor-keyword-color);
  font-weight: bold;
}
.cm-s-jupyter span.cm-atom {
  color: var(--jp-mirror-editor-atom-color);
}
.cm-s-jupyter span.cm-number {
  color: var(--jp-mirror-editor-number-color);
}
.cm-s-jupyter span.cm-def {
  color: var(--jp-mirror-editor-def-color);
}
.cm-s-jupyter span.cm-variable {
  color: var(--jp-mirror-editor-variable-color);
}
.cm-s-jupyter span.cm-variable-2 {
  color: var(--jp-mirror-editor-variable-2-color);
}
.cm-s-jupyter span.cm-variable-3 {
  color: var(--jp-mirror-editor-variable-3-color);
}
.cm-s-jupyter span.cm-punctuation {
  color: var(--jp-mirror-editor-punctuation-color);
}
.cm-s-jupyter span.cm-property {
  color: var(--jp-mirror-editor-property-color);
}
.cm-s-jupyter span.cm-operator {
  color: var(--jp-mirror-editor-operator-color);
  font-weight: bold;
}
.cm-s-jupyter span.cm-comment {
  color: var(--jp-mirror-editor-comment-color);
  font-style: italic;
}
.cm-s-jupyter span.cm-string {
  color: var(--jp-mirror-editor-string-color);
}
.cm-s-jupyter span.cm-string-2 {
  color: var(--jp-mirror-editor-string-2-color);
}
.cm-s-jupyter span.cm-meta {
  color: var(--jp-mirror-editor-meta-color);
}
.cm-s-jupyter span.cm-qualifier {
  color: var(--jp-mirror-editor-qualifier-color);
}
.cm-s-jupyter span.cm-builtin {
  color: var(--jp-mirror-editor-builtin-color);
}
.cm-s-jupyter span.cm-bracket {
  color: var(--jp-mirror-editor-bracket-color);
}
.cm-s-jupyter span.cm-tag {
  color: var(--jp-mirror-editor-tag-color);
}
.cm-s-jupyter span.cm-attribute {
  color: var(--jp-mirror-editor-attribute-color);
}
.cm-s-jupyter span.cm-header {
  color: var(--jp-mirror-editor-header-color);
}
.cm-s-jupyter span.cm-quote {
  color: var(--jp-mirror-editor-quote-color);
}
.cm-s-jupyter span.cm-link {
  color: var(--jp-mirror-editor-link-color);
}
.cm-s-jupyter span.cm-error {
  color: var(--jp-mirror-editor-error-color);
}
.cm-s-jupyter span.cm-hr {
  color: #999;
}

.cm-s-jupyter span.cm-tab {
  background: url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADAAAAAMCAYAAAAkuj5RAAAAAXNSR0IArs4c6QAAAGFJREFUSMft1LsRQFAQheHPowAKoACx3IgEKtaEHujDjORSgWTH/ZOdnZOcM/sgk/kFFWY0qV8foQwS4MKBCS3qR6ixBJvElOobYAtivseIE120FaowJPN75GMu8j/LfMwNjh4HUpwg4LUAAAAASUVORK5CYII=);
  background-position: right;
  background-repeat: no-repeat;
}

.cm-s-jupyter .CodeMirror-activeline-background,
.cm-s-jupyter .CodeMirror-gutter {
  background-color: var(--jp-layout-color2);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| RenderedText
|----------------------------------------------------------------------------*/

.jp-RenderedText {
  text-align: left;
  padding-left: var(--jp-code-padding);
  line-height: var(--jp-code-line-height);
  font-family: var(--jp-code-font-family);
}

.jp-RenderedText pre,
.jp-RenderedJavaScript pre,
.jp-RenderedHTMLCommon pre {
  color: var(--jp-content-font-color1);
  font-size: var(--jp-code-font-size);
  border: none;
  margin: 0px;
  padding: 0px;
  line-height: normal;
}

.jp-RenderedText pre a:link {
  text-decoration: none;
  color: var(--jp-content-link-color);
}
.jp-RenderedText pre a:hover {
  text-decoration: underline;
  color: var(--jp-content-link-color);
}
.jp-RenderedText pre a:visited {
  text-decoration: none;
  color: var(--jp-content-link-color);
}

/* console foregrounds and backgrounds */
.jp-RenderedText pre .ansi-black-fg {
  color: #3e424d;
}
.jp-RenderedText pre .ansi-red-fg {
  color: #e75c58;
}
.jp-RenderedText pre .ansi-green-fg {
  color: #00a250;
}
.jp-RenderedText pre .ansi-yellow-fg {
  color: #ddb62b;
}
.jp-RenderedText pre .ansi-blue-fg {
  color: #208ffb;
}
.jp-RenderedText pre .ansi-magenta-fg {
  color: #d160c4;
}
.jp-RenderedText pre .ansi-cyan-fg {
  color: #60c6c8;
}
.jp-RenderedText pre .ansi-white-fg {
  color: #c5c1b4;
}

.jp-RenderedText pre .ansi-black-bg {
  background-color: #3e424d;
}
.jp-RenderedText pre .ansi-red-bg {
  background-color: #e75c58;
}
.jp-RenderedText pre .ansi-green-bg {
  background-color: #00a250;
}
.jp-RenderedText pre .ansi-yellow-bg {
  background-color: #ddb62b;
}
.jp-RenderedText pre .ansi-blue-bg {
  background-color: #208ffb;
}
.jp-RenderedText pre .ansi-magenta-bg {
  background-color: #d160c4;
}
.jp-RenderedText pre .ansi-cyan-bg {
  background-color: #60c6c8;
}
.jp-RenderedText pre .ansi-white-bg {
  background-color: #c5c1b4;
}

.jp-RenderedText pre .ansi-black-intense-fg {
  color: #282c36;
}
.jp-RenderedText pre .ansi-red-intense-fg {
  color: #b22b31;
}
.jp-RenderedText pre .ansi-green-intense-fg {
  color: #007427;
}
.jp-RenderedText pre .ansi-yellow-intense-fg {
  color: #b27d12;
}
.jp-RenderedText pre .ansi-blue-intense-fg {
  color: #0065ca;
}
.jp-RenderedText pre .ansi-magenta-intense-fg {
  color: #a03196;
}
.jp-RenderedText pre .ansi-cyan-intense-fg {
  color: #258f8f;
}
.jp-RenderedText pre .ansi-white-intense-fg {
  color: #a1a6b2;
}

.jp-RenderedText pre .ansi-black-intense-bg {
  background-color: #282c36;
}
.jp-RenderedText pre .ansi-red-intense-bg {
  background-color: #b22b31;
}
.jp-RenderedText pre .ansi-green-intense-bg {
  background-color: #007427;
}
.jp-RenderedText pre .ansi-yellow-intense-bg {
  background-color: #b27d12;
}
.jp-RenderedText pre .ansi-blue-intense-bg {
  background-color: #0065ca;
}
.jp-RenderedText pre .ansi-magenta-intense-bg {
  background-color: #a03196;
}
.jp-RenderedText pre .ansi-cyan-intense-bg {
  background-color: #258f8f;
}
.jp-RenderedText pre .ansi-white-intense-bg {
  background-color: #a1a6b2;
}

.jp-RenderedText pre .ansi-default-inverse-fg {
  color: var(--jp-ui-inverse-font-color0);
}
.jp-RenderedText pre .ansi-default-inverse-bg {
  background-color: var(--jp-inverse-layout-color0);
}

.jp-RenderedText pre .ansi-bold {
  font-weight: bold;
}
.jp-RenderedText pre .ansi-underline {
  text-decoration: underline;
}

.jp-RenderedText[data-mime-type='application/vnd.jupyter.stderr'] {
  background: var(--jp-rendermime-error-background);
  padding-top: var(--jp-code-padding);
}

/*-----------------------------------------------------------------------------
| RenderedLatex
|----------------------------------------------------------------------------*/

.jp-RenderedLatex {
  color: var(--jp-content-font-color1);
  font-size: var(--jp-content-font-size1);
  line-height: var(--jp-content-line-height);
}

/* Left-justify outputs.*/
.jp-OutputArea-output.jp-RenderedLatex {
  padding: var(--jp-code-padding);
  text-align: left;
}

/*-----------------------------------------------------------------------------
| RenderedHTML
|----------------------------------------------------------------------------*/

.jp-RenderedHTMLCommon {
  color: var(--jp-content-font-color1);
  font-family: var(--jp-content-font-family);
  font-size: var(--jp-content-font-size1);
  line-height: var(--jp-content-line-height);
  /* Give a bit more R padding on Markdown text to keep line lengths reasonable */
  padding-right: 20px;
}

.jp-RenderedHTMLCommon em {
  font-style: italic;
}

.jp-RenderedHTMLCommon strong {
  font-weight: bold;
}

.jp-RenderedHTMLCommon u {
  text-decoration: underline;
}

.jp-RenderedHTMLCommon a:link {
  text-decoration: none;
  color: var(--jp-content-link-color);
}

.jp-RenderedHTMLCommon a:hover {
  text-decoration: underline;
  color: var(--jp-content-link-color);
}

.jp-RenderedHTMLCommon a:visited {
  text-decoration: none;
  color: var(--jp-content-link-color);
}

/* Headings */

.jp-RenderedHTMLCommon h1,
.jp-RenderedHTMLCommon h2,
.jp-RenderedHTMLCommon h3,
.jp-RenderedHTMLCommon h4,
.jp-RenderedHTMLCommon h5,
.jp-RenderedHTMLCommon h6 {
  line-height: var(--jp-content-heading-line-height);
  font-weight: var(--jp-content-heading-font-weight);
  font-style: normal;
  margin: var(--jp-content-heading-margin-top) 0
    var(--jp-content-heading-margin-bottom) 0;
}

.jp-RenderedHTMLCommon h1:first-child,
.jp-RenderedHTMLCommon h2:first-child,
.jp-RenderedHTMLCommon h3:first-child,
.jp-RenderedHTMLCommon h4:first-child,
.jp-RenderedHTMLCommon h5:first-child,
.jp-RenderedHTMLCommon h6:first-child {
  margin-top: calc(0.5 * var(--jp-content-heading-margin-top));
}

.jp-RenderedHTMLCommon h1:last-child,
.jp-RenderedHTMLCommon h2:last-child,
.jp-RenderedHTMLCommon h3:last-child,
.jp-RenderedHTMLCommon h4:last-child,
.jp-RenderedHTMLCommon h5:last-child,
.jp-RenderedHTMLCommon h6:last-child {
  margin-bottom: calc(0.5 * var(--jp-content-heading-margin-bottom));
}

.jp-RenderedHTMLCommon h1 {
  font-size: var(--jp-content-font-size5);
}

.jp-RenderedHTMLCommon h2 {
  font-size: var(--jp-content-font-size4);
}

.jp-RenderedHTMLCommon h3 {
  font-size: var(--jp-content-font-size3);
}

.jp-RenderedHTMLCommon h4 {
  font-size: var(--jp-content-font-size2);
}

.jp-RenderedHTMLCommon h5 {
  font-size: var(--jp-content-font-size1);
}

.jp-RenderedHTMLCommon h6 {
  font-size: var(--jp-content-font-size0);
}

/* Lists */

.jp-RenderedHTMLCommon ul:not(.list-inline),
.jp-RenderedHTMLCommon ol:not(.list-inline) {
  padding-left: 2em;
}

.jp-RenderedHTMLCommon ul {
  list-style: disc;
}

.jp-RenderedHTMLCommon ul ul {
  list-style: square;
}

.jp-RenderedHTMLCommon ul ul ul {
  list-style: circle;
}

.jp-RenderedHTMLCommon ol {
  list-style: decimal;
}

.jp-RenderedHTMLCommon ol ol {
  list-style: upper-alpha;
}

.jp-RenderedHTMLCommon ol ol ol {
  list-style: lower-alpha;
}

.jp-RenderedHTMLCommon ol ol ol ol {
  list-style: lower-roman;
}

.jp-RenderedHTMLCommon ol ol ol ol ol {
  list-style: decimal;
}

.jp-RenderedHTMLCommon ol,
.jp-RenderedHTMLCommon ul {
  margin-bottom: 1em;
}

.jp-RenderedHTMLCommon ul ul,
.jp-RenderedHTMLCommon ul ol,
.jp-RenderedHTMLCommon ol ul,
.jp-RenderedHTMLCommon ol ol {
  margin-bottom: 0em;
}

.jp-RenderedHTMLCommon hr {
  color: var(--jp-border-color2);
  background-color: var(--jp-border-color1);
  margin-top: 1em;
  margin-bottom: 1em;
}

.jp-RenderedHTMLCommon > pre {
  margin: 1.5em 2em;
}

.jp-RenderedHTMLCommon pre,
.jp-RenderedHTMLCommon code {
  border: 0;
  background-color: var(--jp-layout-color0);
  color: var(--jp-content-font-color1);
  font-family: var(--jp-code-font-family);
  font-size: inherit;
  line-height: var(--jp-code-line-height);
  padding: 0;
  white-space: pre-wrap;
}

.jp-RenderedHTMLCommon :not(pre) > code {
  background-color: var(--jp-layout-color2);
  padding: 1px 5px;
}

/* Tables */

.jp-RenderedHTMLCommon table {
  border-collapse: collapse;
  border-spacing: 0;
  border: none;
  color: var(--jp-ui-font-color1);
  font-size: 12px;
  table-layout: fixed;
  margin-left: auto;
  margin-right: auto;
}

.jp-RenderedHTMLCommon thead {
  border-bottom: var(--jp-border-width) solid var(--jp-border-color1);
  vertical-align: bottom;
}

.jp-RenderedHTMLCommon td,
.jp-RenderedHTMLCommon th,
.jp-RenderedHTMLCommon tr {
  vertical-align: middle;
  padding: 0.5em 0.5em;
  line-height: normal;
  white-space: normal;
  max-width: none;
  border: none;
}

.jp-RenderedMarkdown.jp-RenderedHTMLCommon td,
.jp-RenderedMarkdown.jp-RenderedHTMLCommon th {
  max-width: none;
}

:not(.jp-RenderedMarkdown).jp-RenderedHTMLCommon td,
:not(.jp-RenderedMarkdown).jp-RenderedHTMLCommon th,
:not(.jp-RenderedMarkdown).jp-RenderedHTMLCommon tr {
  text-align: right;
}

.jp-RenderedHTMLCommon th {
  font-weight: bold;
}

.jp-RenderedHTMLCommon tbody tr:nth-child(odd) {
  background: var(--jp-layout-color0);
}

.jp-RenderedHTMLCommon tbody tr:nth-child(even) {
  background: var(--jp-rendermime-table-row-background);
}

.jp-RenderedHTMLCommon tbody tr:hover {
  background: var(--jp-rendermime-table-row-hover-background);
}

.jp-RenderedHTMLCommon table {
  margin-bottom: 1em;
}

.jp-RenderedHTMLCommon p {
  text-align: left;
  margin: 0px;
}

.jp-RenderedHTMLCommon p {
  margin-bottom: 1em;
}

.jp-RenderedHTMLCommon img {
  -moz-force-broken-image-icon: 1;
}

/* Restrict to direct children as other images could be nested in other content. */
.jp-RenderedHTMLCommon > img {
  display: block;
  margin-left: 0;
  margin-right: 0;
  margin-bottom: 1em;
}

/* Change color behind transparent images if they need it... */
[data-jp-theme-light='false'] .jp-RenderedImage img.jp-needs-light-background {
  background-color: var(--jp-inverse-layout-color1);
}
[data-jp-theme-light='true'] .jp-RenderedImage img.jp-needs-dark-background {
  background-color: var(--jp-inverse-layout-color1);
}
/* ...or leave it untouched if they don't */
[data-jp-theme-light='false'] .jp-RenderedImage img.jp-needs-dark-background {
}
[data-jp-theme-light='true'] .jp-RenderedImage img.jp-needs-light-background {
}

.jp-RenderedHTMLCommon img,
.jp-RenderedImage img,
.jp-RenderedHTMLCommon svg,
.jp-RenderedSVG svg {
  max-width: 100%;
  height: auto;
}

.jp-RenderedHTMLCommon img.jp-mod-unconfined,
.jp-RenderedImage img.jp-mod-unconfined,
.jp-RenderedHTMLCommon svg.jp-mod-unconfined,
.jp-RenderedSVG svg.jp-mod-unconfined {
  max-width: none;
}

.jp-RenderedHTMLCommon .alert {
  padding: var(--jp-notebook-padding);
  border: var(--jp-border-width) solid transparent;
  border-radius: var(--jp-border-radius);
  margin-bottom: 1em;
}

.jp-RenderedHTMLCommon .alert-info {
  color: var(--jp-info-color0);
  background-color: var(--jp-info-color3);
  border-color: var(--jp-info-color2);
}
.jp-RenderedHTMLCommon .alert-info hr {
  border-color: var(--jp-info-color3);
}
.jp-RenderedHTMLCommon .alert-info > p:last-child,
.jp-RenderedHTMLCommon .alert-info > ul:last-child {
  margin-bottom: 0;
}

.jp-RenderedHTMLCommon .alert-warning {
  color: var(--jp-warn-color0);
  background-color: var(--jp-warn-color3);
  border-color: var(--jp-warn-color2);
}
.jp-RenderedHTMLCommon .alert-warning hr {
  border-color: var(--jp-warn-color3);
}
.jp-RenderedHTMLCommon .alert-warning > p:last-child,
.jp-RenderedHTMLCommon .alert-warning > ul:last-child {
  margin-bottom: 0;
}

.jp-RenderedHTMLCommon .alert-success {
  color: var(--jp-success-color0);
  background-color: var(--jp-success-color3);
  border-color: var(--jp-success-color2);
}
.jp-RenderedHTMLCommon .alert-success hr {
  border-color: var(--jp-success-color3);
}
.jp-RenderedHTMLCommon .alert-success > p:last-child,
.jp-RenderedHTMLCommon .alert-success > ul:last-child {
  margin-bottom: 0;
}

.jp-RenderedHTMLCommon .alert-danger {
  color: var(--jp-error-color0);
  background-color: var(--jp-error-color3);
  border-color: var(--jp-error-color2);
}
.jp-RenderedHTMLCommon .alert-danger hr {
  border-color: var(--jp-error-color3);
}
.jp-RenderedHTMLCommon .alert-danger > p:last-child,
.jp-RenderedHTMLCommon .alert-danger > ul:last-child {
  margin-bottom: 0;
}

.jp-RenderedHTMLCommon blockquote {
  margin: 1em 2em;
  padding: 0 1em;
  border-left: 5px solid var(--jp-border-color2);
}

a.jp-InternalAnchorLink {
  visibility: hidden;
  margin-left: 8px;
  color: var(--md-blue-800);
}

h1:hover .jp-InternalAnchorLink,
h2:hover .jp-InternalAnchorLink,
h3:hover .jp-InternalAnchorLink,
h4:hover .jp-InternalAnchorLink,
h5:hover .jp-InternalAnchorLink,
h6:hover .jp-InternalAnchorLink {
  visibility: visible;
}

.jp-RenderedHTMLCommon kbd {
  background-color: var(--jp-rendermime-table-row-background);
  border: 1px solid var(--jp-border-color0);
  border-bottom-color: var(--jp-border-color2);
  border-radius: 3px;
  box-shadow: inset 0 -1px 0 rgba(0, 0, 0, 0.25);
  display: inline-block;
  font-size: 0.8em;
  line-height: 1em;
  padding: 0.2em 0.5em;
}

/* Most direct children of .jp-RenderedHTMLCommon have a margin-bottom of 1.0.
 * At the bottom of cells this is a bit too much as there is also spacing
 * between cells. Going all the way to 0 gets too tight between markdown and
 * code cells.
 */
.jp-RenderedHTMLCommon > *:last-child {
  margin-bottom: 0.5em;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-MimeDocument {
  outline: none;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Variables
|----------------------------------------------------------------------------*/

:root {
  --jp-private-filebrowser-button-height: 28px;
  --jp-private-filebrowser-button-width: 48px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-FileBrowser {
  display: flex;
  flex-direction: column;
  color: var(--jp-ui-font-color1);
  background: var(--jp-layout-color1);
  /* This is needed so that all font sizing of children done in ems is
   * relative to this base size */
  font-size: var(--jp-ui-font-size1);
}

.jp-FileBrowser-toolbar.jp-Toolbar {
  border-bottom: none;
  height: auto;
  margin: var(--jp-toolbar-header-margin);
  box-shadow: none;
}

.jp-BreadCrumbs {
  flex: 0 0 auto;
  margin: 4px 12px;
}

.jp-BreadCrumbs-item {
  margin: 0px 2px;
  padding: 0px 2px;
  border-radius: var(--jp-border-radius);
  cursor: pointer;
}

.jp-BreadCrumbs-item:hover {
  background-color: var(--jp-layout-color2);
}

.jp-BreadCrumbs-item:first-child {
  margin-left: 0px;
}

.jp-BreadCrumbs-item.jp-mod-dropTarget {
  background-color: var(--jp-brand-color2);
  opacity: 0.7;
}

/*-----------------------------------------------------------------------------
| Buttons
|----------------------------------------------------------------------------*/

.jp-FileBrowser-toolbar.jp-Toolbar {
  padding: 0px;
}

.jp-FileBrowser-toolbar.jp-Toolbar {
  justify-content: space-evenly;
}

.jp-FileBrowser-toolbar.jp-Toolbar .jp-Toolbar-item {
  flex: 1;
}

.jp-FileBrowser-toolbar.jp-Toolbar .jp-ToolbarButtonComponent {
  width: 100%;
}

/*-----------------------------------------------------------------------------
| DirListing
|----------------------------------------------------------------------------*/

.jp-DirListing {
  flex: 1 1 auto;
  display: flex;
  flex-direction: column;
  outline: 0;
}

.jp-DirListing-header {
  flex: 0 0 auto;
  display: flex;
  flex-direction: row;
  overflow: hidden;
  border-top: var(--jp-border-width) solid var(--jp-border-color2);
  border-bottom: var(--jp-border-width) solid var(--jp-border-color1);
  box-shadow: var(--jp-toolbar-box-shadow);
  z-index: 2;
}

.jp-DirListing-headerItem {
  padding: 4px 12px 2px 12px;
  font-weight: 500;
}

.jp-DirListing-headerItem:hover {
  background: var(--jp-layout-color2);
}

.jp-DirListing-headerItem.jp-id-name {
  flex: 1 0 84px;
}

.jp-DirListing-headerItem.jp-id-modified {
  flex: 0 0 112px;
  border-left: var(--jp-border-width) solid var(--jp-border-color2);
  text-align: right;
}

.jp-DirListing-narrow .jp-id-modified,
.jp-DirListing-narrow .jp-DirListing-itemModified {
  display: none;
}

.jp-DirListing-headerItem.jp-mod-selected {
  font-weight: 600;
}

/* increase specificity to override bundled default */
.jp-DirListing-content {
  flex: 1 1 auto;
  margin: 0;
  padding: 0;
  list-style-type: none;
  overflow: auto;
  background-color: var(--jp-layout-color1);
}

/* Style the directory listing content when a user drops a file to upload */
.jp-DirListing.jp-mod-native-drop .jp-DirListing-content {
  outline: 5px dashed rgba(128, 128, 128, 0.5);
  outline-offset: -10px;
  cursor: copy;
}

.jp-DirListing-item {
  display: flex;
  flex-direction: row;
  padding: 4px 12px;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}

.jp-DirListing-item.jp-mod-selected {
  color: white;
  background: var(--jp-brand-color1);
}

.jp-DirListing-item.jp-mod-dropTarget {
  background: var(--jp-brand-color3);
}

.jp-DirListing-item:hover:not(.jp-mod-selected) {
  background: var(--jp-layout-color2);
}

.jp-DirListing-itemIcon {
  flex: 0 0 20px;
  margin-right: 4px;
}

.jp-DirListing-itemText {
  flex: 1 0 64px;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  user-select: none;
}

.jp-DirListing-itemModified {
  flex: 0 0 125px;
  text-align: right;
}

.jp-DirListing-editor {
  flex: 1 0 64px;
  outline: none;
  border: none;
}

.jp-DirListing-item.jp-mod-running .jp-DirListing-itemIcon:before {
  color: limegreen;
  content: '\25CF';
  font-size: 8px;
  position: absolute;
  left: -8px;
}

.jp-DirListing-item.lm-mod-drag-image,
.jp-DirListing-item.jp-mod-selected.lm-mod-drag-image {
  font-size: var(--jp-ui-font-size1);
  padding-left: 4px;
  margin-left: 4px;
  width: 160px;
  background-color: var(--jp-ui-inverse-font-color2);
  box-shadow: var(--jp-elevation-z2);
  border-radius: 0px;
  color: var(--jp-ui-font-color1);
  transform: translateX(-40%) translateY(-58%);
}

.jp-DirListing-deadSpace {
  flex: 1 1 auto;
  margin: 0;
  padding: 0;
  list-style-type: none;
  overflow: auto;
  background-color: var(--jp-layout-color1);
}

.jp-Document {
  min-width: 120px;
  min-height: 120px;
  outline: none;
}

.jp-FileDialog.jp-mod-conflict input {
  color: red;
}

.jp-FileDialog .jp-new-name-title {
  margin-top: 12px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Private CSS variables
|----------------------------------------------------------------------------*/

:root {
}

/*-----------------------------------------------------------------------------
| Main OutputArea
| OutputArea has a list of Outputs
|----------------------------------------------------------------------------*/

.jp-OutputArea {
  overflow-y: auto;
}

.jp-OutputArea-child {
  display: flex;
  flex-direction: row;
}

.jp-OutputPrompt {
  flex: 0 0 var(--jp-cell-prompt-width);
  color: var(--jp-cell-outprompt-font-color);
  font-family: var(--jp-cell-prompt-font-family);
  padding: var(--jp-code-padding);
  letter-spacing: var(--jp-cell-prompt-letter-spacing);
  line-height: var(--jp-code-line-height);
  font-size: var(--jp-code-font-size);
  border: var(--jp-border-width) solid transparent;
  opacity: var(--jp-cell-prompt-opacity);
  /* Right align prompt text, don't wrap to handle large prompt numbers */
  text-align: right;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  /* Disable text selection */
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}

.jp-OutputArea-output {
  height: auto;
  overflow: auto;
  user-select: text;
  -moz-user-select: text;
  -webkit-user-select: text;
  -ms-user-select: text;
}

.jp-OutputArea-child .jp-OutputArea-output {
  flex-grow: 1;
  flex-shrink: 1;
}

/**
 * Isolated output.
 */
.jp-OutputArea-output.jp-mod-isolated {
  width: 100%;
  display: block;
}

/*
When drag events occur, `p-mod-override-cursor` is added to the body.
Because iframes steal all cursor events, the following two rules are necessary
to suppress pointer events while resize drags are occurring. There may be a
better solution to this problem.
*/
body.lm-mod-override-cursor .jp-OutputArea-output.jp-mod-isolated {
  position: relative;
}

body.lm-mod-override-cursor .jp-OutputArea-output.jp-mod-isolated:before {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: transparent;
}

/* pre */

.jp-OutputArea-output pre {
  border: none;
  margin: 0px;
  padding: 0px;
  overflow-x: auto;
  overflow-y: auto;
  word-break: break-all;
  word-wrap: break-word;
  white-space: pre-wrap;
}

/* tables */

.jp-OutputArea-output.jp-RenderedHTMLCommon table {
  margin-left: 0;
  margin-right: 0;
}

/* description lists */

.jp-OutputArea-output dl,
.jp-OutputArea-output dt,
.jp-OutputArea-output dd {
  display: block;
}

.jp-OutputArea-output dl {
  width: 100%;
  overflow: hidden;
  padding: 0;
  margin: 0;
}

.jp-OutputArea-output dt {
  font-weight: bold;
  float: left;
  width: 20%;
  padding: 0;
  margin: 0;
}

.jp-OutputArea-output dd {
  float: left;
  width: 80%;
  padding: 0;
  margin: 0;
}

/* Hide the gutter in case of
 *  - nested output areas (e.g. in the case of output widgets)
 *  - mirrored output areas
 */
.jp-OutputArea .jp-OutputArea .jp-OutputArea-prompt {
  display: none;
}

/*-----------------------------------------------------------------------------
| executeResult is added to any Output-result for the display of the object
| returned by a cell
|----------------------------------------------------------------------------*/

.jp-OutputArea-output.jp-OutputArea-executeResult {
  margin-left: 0px;
  flex: 1 1 auto;
}

.jp-OutputArea-executeResult.jp-RenderedText {
  padding-top: var(--jp-code-padding);
}

/*-----------------------------------------------------------------------------
| The Stdin output
|----------------------------------------------------------------------------*/

.jp-OutputArea-stdin {
  line-height: var(--jp-code-line-height);
  padding-top: var(--jp-code-padding);
  display: flex;
}

.jp-Stdin-prompt {
  color: var(--jp-content-font-color0);
  padding-right: var(--jp-code-padding);
  vertical-align: baseline;
  flex: 0 0 auto;
}

.jp-Stdin-input {
  font-family: var(--jp-code-font-family);
  font-size: inherit;
  color: inherit;
  background-color: inherit;
  width: 42%;
  min-width: 200px;
  /* make sure input baseline aligns with prompt */
  vertical-align: baseline;
  /* padding + margin = 0.5em between prompt and cursor */
  padding: 0em 0.25em;
  margin: 0em 0.25em;
  flex: 0 0 70%;
}

.jp-Stdin-input:focus {
  box-shadow: none;
}

/*-----------------------------------------------------------------------------
| Output Area View
|----------------------------------------------------------------------------*/

.jp-LinkedOutputView .jp-OutputArea {
  height: 100%;
  display: block;
}

.jp-LinkedOutputView .jp-OutputArea-output:only-child {
  height: 100%;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-Collapser {
  flex: 0 0 var(--jp-cell-collapser-width);
  padding: 0px;
  margin: 0px;
  border: none;
  outline: none;
  background: transparent;
  border-radius: var(--jp-border-radius);
  opacity: 1;
}

.jp-Collapser-child {
  display: block;
  width: 100%;
  box-sizing: border-box;
  /* height: 100% doesn't work because the height of its parent is computed from content */
  position: absolute;
  top: 0px;
  bottom: 0px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Header/Footer
|----------------------------------------------------------------------------*/

/* Hidden by zero height by default */
.jp-CellHeader,
.jp-CellFooter {
  height: 0px;
  width: 100%;
  padding: 0px;
  margin: 0px;
  border: none;
  outline: none;
  background: transparent;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Input
|----------------------------------------------------------------------------*/

/* All input areas */
.jp-InputArea {
  display: flex;
  flex-direction: row;
}

.jp-InputArea-editor {
  flex: 1 1 auto;
}

.jp-InputArea-editor {
  /* This is the non-active, default styling */
  border: var(--jp-border-width) solid var(--jp-cell-editor-border-color);
  border-radius: 0px;
  background: var(--jp-cell-editor-background);
}

.jp-InputPrompt {
  flex: 0 0 var(--jp-cell-prompt-width);
  color: var(--jp-cell-inprompt-font-color);
  font-family: var(--jp-cell-prompt-font-family);
  padding: var(--jp-code-padding);
  letter-spacing: var(--jp-cell-prompt-letter-spacing);
  opacity: var(--jp-cell-prompt-opacity);
  line-height: var(--jp-code-line-height);
  font-size: var(--jp-code-font-size);
  border: var(--jp-border-width) solid transparent;
  opacity: var(--jp-cell-prompt-opacity);
  /* Right align prompt text, don't wrap to handle large prompt numbers */
  text-align: right;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  /* Disable text selection */
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Placeholder
|----------------------------------------------------------------------------*/

.jp-Placeholder {
  display: flex;
  flex-direction: row;
  flex: 1 1 auto;
}

.jp-Placeholder-prompt {
  box-sizing: border-box;
}

.jp-Placeholder-content {
  flex: 1 1 auto;
  border: none;
  background: transparent;
  height: 20px;
  box-sizing: border-box;
}

.jp-Placeholder-content .jp-MoreHorizIcon {
  width: 32px;
  height: 16px;
  border: 1px solid transparent;
  border-radius: var(--jp-border-radius);
}

.jp-Placeholder-content .jp-MoreHorizIcon:hover {
  border: 1px solid var(--jp-border-color1);
  box-shadow: 0px 0px 2px 0px rgba(0, 0, 0, 0.25);
  background-color: var(--jp-layout-color0);
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Private CSS variables
|----------------------------------------------------------------------------*/

:root {
  --jp-private-cell-scrolling-output-offset: 5px;
}

/*-----------------------------------------------------------------------------
| Cell
|----------------------------------------------------------------------------*/

.jp-Cell {
  padding: var(--jp-cell-padding);
  margin: 0px;
  border: none;
  outline: none;
  background: transparent;
}

/*-----------------------------------------------------------------------------
| Common input/output
|----------------------------------------------------------------------------*/

.jp-Cell-inputWrapper,
.jp-Cell-outputWrapper {
  display: flex;
  flex-direction: row;
  padding: 0px;
  margin: 0px;
  /* Added to reveal the box-shadow on the input and output collapsers. */
  overflow: visible;
}

/* Only input/output areas inside cells */
.jp-Cell-inputArea,
.jp-Cell-outputArea {
  flex: 1 1 auto;
}

/*-----------------------------------------------------------------------------
| Collapser
|----------------------------------------------------------------------------*/

/* Make the output collapser disappear when there is not output, but do so
 * in a manner that leaves it in the layout and preserves its width.
 */
.jp-Cell.jp-mod-noOutputs .jp-Cell-outputCollapser {
  border: none !important;
  background: transparent !important;
}

.jp-Cell:not(.jp-mod-noOutputs) .jp-Cell-outputCollapser {
  min-height: var(--jp-cell-collapser-min-height);
}

/*-----------------------------------------------------------------------------
| Output
|----------------------------------------------------------------------------*/

/* Put a space between input and output when there IS output */
.jp-Cell:not(.jp-mod-noOutputs) .jp-Cell-outputWrapper {
  margin-top: 5px;
}

/* Text output with the Out[] prompt needs a top padding to match the
 * alignment of the Out[] prompt itself.
 */
.jp-OutputArea-executeResult .jp-RenderedText.jp-OutputArea-output {
  padding-top: var(--jp-code-padding);
}

.jp-CodeCell.jp-mod-outputsScrolled .jp-Cell-outputArea {
  overflow-y: auto;
  max-height: 200px;
  box-shadow: inset 0 0 6px 2px rgba(0, 0, 0, 0.3);
  margin-left: var(--jp-private-cell-scrolling-output-offset);
}

.jp-CodeCell.jp-mod-outputsScrolled .jp-OutputArea-prompt {
  flex: 0 0
    calc(
      var(--jp-cell-prompt-width) -
        var(--jp-private-cell-scrolling-output-offset)
    );
}

/*-----------------------------------------------------------------------------
| CodeCell
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| MarkdownCell
|----------------------------------------------------------------------------*/

.jp-MarkdownOutput {
  flex: 1 1 auto;
  margin-top: 0;
  margin-bottom: 0;
  padding-left: var(--jp-code-padding);
}

.jp-MarkdownOutput.jp-RenderedHTMLCommon {
  overflow: auto;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Variables
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------

/*-----------------------------------------------------------------------------
| Styles
|----------------------------------------------------------------------------*/

.jp-NotebookPanel-toolbar {
  padding: 2px;
}

.jp-Toolbar-item.jp-Notebook-toolbarCellType .jp-select-wrapper.jp-mod-focused {
  border: none;
  box-shadow: none;
}

.jp-Notebook-toolbarCellTypeDropdown select {
  height: 24px;
  font-size: var(--jp-ui-font-size1);
  line-height: 14px;
  border-radius: 0;
  display: block;
}

.jp-Notebook-toolbarCellTypeDropdown span {
  top: 5px !important;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Private CSS variables
|----------------------------------------------------------------------------*/

:root {
  --jp-private-notebook-dragImage-width: 304px;
  --jp-private-notebook-dragImage-height: 36px;
  --jp-private-notebook-selected-color: var(--md-blue-400);
  --jp-private-notebook-active-color: var(--md-green-400);
}

/*-----------------------------------------------------------------------------
| Imports
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Notebook
|----------------------------------------------------------------------------*/

.jp-NotebookPanel {
  display: block;
  height: 100%;
}

.jp-NotebookPanel.jp-Document {
  min-width: 240px;
  min-height: 120px;
}

.jp-Notebook {
  padding: var(--jp-notebook-padding);
  outline: none;
  overflow: auto;
  background: var(--jp-layout-color0);
}

.jp-Notebook.jp-mod-scrollPastEnd::after {
  display: block;
  content: '';
  min-height: var(--jp-notebook-scroll-padding);
}

.jp-Notebook .jp-Cell {
  overflow: visible;
}

.jp-Notebook .jp-Cell .jp-InputPrompt {
  cursor: move;
}

/*-----------------------------------------------------------------------------
| Notebook state related styling
|
| The notebook and cells each have states, here are the possibilities:
|
| - Notebook
|   - Command
|   - Edit
| - Cell
|   - None
|   - Active (only one can be active)
|   - Selected (the cells actions are applied to)
|   - Multiselected (when multiple selected, the cursor)
|   - No outputs
|----------------------------------------------------------------------------*/

/* Command or edit modes */

.jp-Notebook .jp-Cell:not(.jp-mod-active) .jp-InputPrompt {
  opacity: var(--jp-cell-prompt-not-active-opacity);
  color: var(--jp-cell-prompt-not-active-font-color);
}

.jp-Notebook .jp-Cell:not(.jp-mod-active) .jp-OutputPrompt {
  opacity: var(--jp-cell-prompt-not-active-opacity);
  color: var(--jp-cell-prompt-not-active-font-color);
}

/* cell is active */
.jp-Notebook .jp-Cell.jp-mod-active .jp-Collapser {
  background: var(--jp-brand-color1);
}

/* collapser is hovered */
.jp-Notebook .jp-Cell .jp-Collapser:hover {
  box-shadow: var(--jp-elevation-z2);
  background: var(--jp-brand-color1);
  opacity: var(--jp-cell-collapser-not-active-hover-opacity);
}

/* cell is active and collapser is hovered */
.jp-Notebook .jp-Cell.jp-mod-active .jp-Collapser:hover {
  background: var(--jp-brand-color0);
  opacity: 1;
}

/* Command mode */

.jp-Notebook.jp-mod-commandMode .jp-Cell.jp-mod-selected {
  background: var(--jp-notebook-multiselected-color);
}

.jp-Notebook.jp-mod-commandMode
  .jp-Cell.jp-mod-active.jp-mod-selected:not(.jp-mod-multiSelected) {
  background: transparent;
}

/* Edit mode */

.jp-Notebook.jp-mod-editMode .jp-Cell.jp-mod-active .jp-InputArea-editor {
  border: var(--jp-border-width) solid var(--jp-cell-editor-active-border-color);
  box-shadow: var(--jp-input-box-shadow);
  background-color: var(--jp-cell-editor-active-background);
}

/*-----------------------------------------------------------------------------
| Notebook drag and drop
|----------------------------------------------------------------------------*/

.jp-Notebook-cell.jp-mod-dropSource {
  opacity: 0.5;
}

.jp-Notebook-cell.jp-mod-dropTarget,
.jp-Notebook.jp-mod-commandMode
  .jp-Notebook-cell.jp-mod-active.jp-mod-selected.jp-mod-dropTarget {
  border-top-color: var(--jp-private-notebook-selected-color);
  border-top-style: solid;
  border-top-width: 2px;
}

.jp-dragImage {
  display: flex;
  flex-direction: row;
  width: var(--jp-private-notebook-dragImage-width);
  height: var(--jp-private-notebook-dragImage-height);
  border: var(--jp-border-width) solid var(--jp-cell-editor-border-color);
  background: var(--jp-cell-editor-background);
  overflow: visible;
}

.jp-dragImage-singlePrompt {
  box-shadow: 2px 2px 4px 0px rgba(0, 0, 0, 0.12);
}

.jp-dragImage .jp-dragImage-content {
  flex: 1 1 auto;
  z-index: 2;
  font-size: var(--jp-code-font-size);
  font-family: var(--jp-code-font-family);
  line-height: var(--jp-code-line-height);
  padding: var(--jp-code-padding);
  border: var(--jp-border-width) solid var(--jp-cell-editor-border-color);
  background: var(--jp-cell-editor-background-color);
  color: var(--jp-content-font-color3);
  text-align: left;
  margin: 4px 4px 4px 0px;
}

.jp-dragImage .jp-dragImage-prompt {
  flex: 0 0 auto;
  min-width: 36px;
  color: var(--jp-cell-inprompt-font-color);
  padding: var(--jp-code-padding);
  padding-left: 12px;
  font-family: var(--jp-cell-prompt-font-family);
  letter-spacing: var(--jp-cell-prompt-letter-spacing);
  line-height: 1.9;
  font-size: var(--jp-code-font-size);
  border: var(--jp-border-width) solid transparent;
}

.jp-dragImage-multipleBack {
  z-index: -1;
  position: absolute;
  height: 32px;
  width: 300px;
  top: 8px;
  left: 8px;
  background: var(--jp-layout-color2);
  border: var(--jp-border-width) solid var(--jp-input-border-color);
  box-shadow: 2px 2px 4px 0px rgba(0, 0, 0, 0.12);
}

/*-----------------------------------------------------------------------------
| Cell toolbar
|----------------------------------------------------------------------------*/

.jp-NotebookTools {
  display: block;
  min-width: var(--jp-sidebar-min-width);
  color: var(--jp-ui-font-color1);
  background: var(--jp-layout-color1);
  /* This is needed so that all font sizing of children done in ems is
    * relative to this base size */
  font-size: var(--jp-ui-font-size1);
  overflow: auto;
}

.jp-NotebookTools-tool {
  padding: 0px 12px 0 12px;
}

.jp-ActiveCellTool {
  padding: 12px;
  background-color: var(--jp-layout-color1);
  border-top: none !important;
}

.jp-ActiveCellTool .jp-InputArea-prompt {
  flex: 0 0 auto;
  padding-left: 0px;
}

.jp-ActiveCellTool .jp-InputArea-editor {
  flex: 1 1 auto;
  background: var(--jp-cell-editor-background);
  border-color: var(--jp-cell-editor-border-color);
}

.jp-ActiveCellTool .jp-InputArea-editor .CodeMirror {
  background: transparent;
}

.jp-MetadataEditorTool {
  flex-direction: column;
  padding: 12px 0px 12px 0px;
}

.jp-RankedPanel > :not(:first-child) {
  margin-top: 12px;
}

.jp-KeySelector select.jp-mod-styled {
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color0);
  border: var(--jp-border-width) solid var(--jp-border-color1);
}

.jp-KeySelector label,
.jp-MetadataEditorTool label {
  line-height: 1.4;
}

/*-----------------------------------------------------------------------------
| Presentation Mode (.jp-mod-presentationMode)
|----------------------------------------------------------------------------*/

.jp-mod-presentationMode .jp-Notebook {
  --jp-content-font-size1: var(--jp-content-presentation-font-size1);
  --jp-code-font-size: var(--jp-code-presentation-font-size);
}

.jp-mod-presentationMode .jp-Notebook .jp-Cell .jp-InputPrompt,
.jp-mod-presentationMode .jp-Notebook .jp-Cell .jp-OutputPrompt {
  flex: 0 0 110px;
}

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/* This file was auto-generated by ensurePackage() in @jupyterlab/buildutils */

/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

</style>

    <style type="text/css">
/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/*
The following CSS variables define the main, public API for styling JupyterLab.
These variables should be used by all plugins wherever possible. In other
words, plugins should not define custom colors, sizes, etc unless absolutely
necessary. This enables users to change the visual theme of JupyterLab
by changing these variables.

Many variables appear in an ordered sequence (0,1,2,3). These sequences
are designed to work well together, so for example, `--jp-border-color1` should
be used with `--jp-layout-color1`. The numbers have the following meanings:

* 0: super-primary, reserved for special emphasis
* 1: primary, most important under normal situations
* 2: secondary, next most important under normal situations
* 3: tertiary, next most important under normal situations

Throughout JupyterLab, we are mostly following principles from Google's
Material Design when selecting colors. We are not, however, following
all of MD as it is not optimized for dense, information rich UIs.
*/

:root {
  /* Elevation
   *
   * We style box-shadows using Material Design's idea of elevation. These particular numbers are taken from here:
   *
   * https://github.com/material-components/material-components-web
   * https://material-components-web.appspot.com/elevation.html
   */

  --jp-shadow-base-lightness: 0;
  --jp-shadow-umbra-color: rgba(
    var(--jp-shadow-base-lightness),
    var(--jp-shadow-base-lightness),
    var(--jp-shadow-base-lightness),
    0.2
  );
  --jp-shadow-penumbra-color: rgba(
    var(--jp-shadow-base-lightness),
    var(--jp-shadow-base-lightness),
    var(--jp-shadow-base-lightness),
    0.14
  );
  --jp-shadow-ambient-color: rgba(
    var(--jp-shadow-base-lightness),
    var(--jp-shadow-base-lightness),
    var(--jp-shadow-base-lightness),
    0.12
  );
  --jp-elevation-z0: none;
  --jp-elevation-z1: 0px 2px 1px -1px var(--jp-shadow-umbra-color),
    0px 1px 1px 0px var(--jp-shadow-penumbra-color),
    0px 1px 3px 0px var(--jp-shadow-ambient-color);
  --jp-elevation-z2: 0px 3px 1px -2px var(--jp-shadow-umbra-color),
    0px 2px 2px 0px var(--jp-shadow-penumbra-color),
    0px 1px 5px 0px var(--jp-shadow-ambient-color);
  --jp-elevation-z4: 0px 2px 4px -1px var(--jp-shadow-umbra-color),
    0px 4px 5px 0px var(--jp-shadow-penumbra-color),
    0px 1px 10px 0px var(--jp-shadow-ambient-color);
  --jp-elevation-z6: 0px 3px 5px -1px var(--jp-shadow-umbra-color),
    0px 6px 10px 0px var(--jp-shadow-penumbra-color),
    0px 1px 18px 0px var(--jp-shadow-ambient-color);
  --jp-elevation-z8: 0px 5px 5px -3px var(--jp-shadow-umbra-color),
    0px 8px 10px 1px var(--jp-shadow-penumbra-color),
    0px 3px 14px 2px var(--jp-shadow-ambient-color);
  --jp-elevation-z12: 0px 7px 8px -4px var(--jp-shadow-umbra-color),
    0px 12px 17px 2px var(--jp-shadow-penumbra-color),
    0px 5px 22px 4px var(--jp-shadow-ambient-color);
  --jp-elevation-z16: 0px 8px 10px -5px var(--jp-shadow-umbra-color),
    0px 16px 24px 2px var(--jp-shadow-penumbra-color),
    0px 6px 30px 5px var(--jp-shadow-ambient-color);
  --jp-elevation-z20: 0px 10px 13px -6px var(--jp-shadow-umbra-color),
    0px 20px 31px 3px var(--jp-shadow-penumbra-color),
    0px 8px 38px 7px var(--jp-shadow-ambient-color);
  --jp-elevation-z24: 0px 11px 15px -7px var(--jp-shadow-umbra-color),
    0px 24px 38px 3px var(--jp-shadow-penumbra-color),
    0px 9px 46px 8px var(--jp-shadow-ambient-color);

  /* Borders
   *
   * The following variables, specify the visual styling of borders in JupyterLab.
   */

  --jp-border-width: 1px;
  --jp-border-color0: var(--md-grey-400);
  --jp-border-color1: var(--md-grey-400);
  --jp-border-color2: var(--md-grey-300);
  --jp-border-color3: var(--md-grey-200);
  --jp-border-radius: 2px;

  /* UI Fonts
   *
   * The UI font CSS variables are used for the typography all of the JupyterLab
   * user interface elements that are not directly user generated content.
   *
   * The font sizing here is done assuming that the body font size of --jp-ui-font-size1
   * is applied to a parent element. When children elements, such as headings, are sized
   * in em all things will be computed relative to that body size.
   */

  --jp-ui-font-scale-factor: 1.2;
  --jp-ui-font-size0: 0.83333em;
  --jp-ui-font-size1: 13px; /* Base font size */
  --jp-ui-font-size2: 1.2em;
  --jp-ui-font-size3: 1.44em;

  --jp-ui-font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica,
    Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol';

  /*
   * Use these font colors against the corresponding main layout colors.
   * In a light theme, these go from dark to light.
   */

  /* Defaults use Material Design specification */
  --jp-ui-font-color0: rgba(0, 0, 0, 1);
  --jp-ui-font-color1: rgba(0, 0, 0, 0.87);
  --jp-ui-font-color2: rgba(0, 0, 0, 0.54);
  --jp-ui-font-color3: rgba(0, 0, 0, 0.38);

  /*
   * Use these against the brand/accent/warn/error colors.
   * These will typically go from light to darker, in both a dark and light theme.
   */

  --jp-ui-inverse-font-color0: rgba(255, 255, 255, 1);
  --jp-ui-inverse-font-color1: rgba(255, 255, 255, 1);
  --jp-ui-inverse-font-color2: rgba(255, 255, 255, 0.7);
  --jp-ui-inverse-font-color3: rgba(255, 255, 255, 0.5);

  /* Content Fonts
   *
   * Content font variables are used for typography of user generated content.
   *
   * The font sizing here is done assuming that the body font size of --jp-content-font-size1
   * is applied to a parent element. When children elements, such as headings, are sized
   * in em all things will be computed relative to that body size.
   */

  --jp-content-line-height: 1.6;
  --jp-content-font-scale-factor: 1.2;
  --jp-content-font-size0: 0.83333em;
  --jp-content-font-size1: 14px; /* Base font size */
  --jp-content-font-size2: 1.2em;
  --jp-content-font-size3: 1.44em;
  --jp-content-font-size4: 1.728em;
  --jp-content-font-size5: 2.0736em;

  /* This gives a magnification of about 125% in presentation mode over normal. */
  --jp-content-presentation-font-size1: 17px;

  --jp-content-heading-line-height: 1;
  --jp-content-heading-margin-top: 1.2em;
  --jp-content-heading-margin-bottom: 0.8em;
  --jp-content-heading-font-weight: 500;

  /* Defaults use Material Design specification */
  --jp-content-font-color0: rgba(0, 0, 0, 1);
  --jp-content-font-color1: rgba(0, 0, 0, 0.87);
  --jp-content-font-color2: rgba(0, 0, 0, 0.54);
  --jp-content-font-color3: rgba(0, 0, 0, 0.38);

  --jp-content-link-color: var(--md-blue-700);

  --jp-content-font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI',
    Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji',
    'Segoe UI Symbol';

  /*
   * Code Fonts
   *
   * Code font variables are used for typography of code and other monospaces content.
   */

  --jp-code-font-size: 13px;
  --jp-code-line-height: 1.3077; /* 17px for 13px base */
  --jp-code-padding: 5px; /* 5px for 13px base, codemirror highlighting needs integer px value */
  --jp-code-font-family-default: Menlo, Consolas, 'DejaVu Sans Mono', monospace;
  --jp-code-font-family: var(--jp-code-font-family-default);

  /* This gives a magnification of about 125% in presentation mode over normal. */
  --jp-code-presentation-font-size: 16px;

  /* may need to tweak cursor width if you change font size */
  --jp-code-cursor-width0: 1.4px;
  --jp-code-cursor-width1: 2px;
  --jp-code-cursor-width2: 4px;

  /* Layout
   *
   * The following are the main layout colors use in JupyterLab. In a light
   * theme these would go from light to dark.
   */

  --jp-layout-color0: white;
  --jp-layout-color1: white;
  --jp-layout-color2: var(--md-grey-200);
  --jp-layout-color3: var(--md-grey-400);
  --jp-layout-color4: var(--md-grey-600);

  /* Inverse Layout
   *
   * The following are the inverse layout colors use in JupyterLab. In a light
   * theme these would go from dark to light.
   */

  --jp-inverse-layout-color0: #111111;
  --jp-inverse-layout-color1: var(--md-grey-900);
  --jp-inverse-layout-color2: var(--md-grey-800);
  --jp-inverse-layout-color3: var(--md-grey-700);
  --jp-inverse-layout-color4: var(--md-grey-600);

  /* Brand/accent */

  --jp-brand-color0: var(--md-blue-700);
  --jp-brand-color1: var(--md-blue-500);
  --jp-brand-color2: var(--md-blue-300);
  --jp-brand-color3: var(--md-blue-100);
  --jp-brand-color4: var(--md-blue-50);

  --jp-accent-color0: var(--md-green-700);
  --jp-accent-color1: var(--md-green-500);
  --jp-accent-color2: var(--md-green-300);
  --jp-accent-color3: var(--md-green-100);

  /* State colors (warn, error, success, info) */

  --jp-warn-color0: var(--md-orange-700);
  --jp-warn-color1: var(--md-orange-500);
  --jp-warn-color2: var(--md-orange-300);
  --jp-warn-color3: var(--md-orange-100);

  --jp-error-color0: var(--md-red-700);
  --jp-error-color1: var(--md-red-500);
  --jp-error-color2: var(--md-red-300);
  --jp-error-color3: var(--md-red-100);

  --jp-success-color0: var(--md-green-700);
  --jp-success-color1: var(--md-green-500);
  --jp-success-color2: var(--md-green-300);
  --jp-success-color3: var(--md-green-100);

  --jp-info-color0: var(--md-cyan-700);
  --jp-info-color1: var(--md-cyan-500);
  --jp-info-color2: var(--md-cyan-300);
  --jp-info-color3: var(--md-cyan-100);

  /* Cell specific styles */

  --jp-cell-padding: 5px;

  --jp-cell-collapser-width: 8px;
  --jp-cell-collapser-min-height: 20px;
  --jp-cell-collapser-not-active-hover-opacity: 0.6;

  --jp-cell-editor-background: var(--md-grey-100);
  --jp-cell-editor-border-color: var(--md-grey-300);
  --jp-cell-editor-box-shadow: inset 0 0 2px var(--md-blue-300);
  --jp-cell-editor-active-background: var(--jp-layout-color0);
  --jp-cell-editor-active-border-color: var(--jp-brand-color1);

  --jp-cell-prompt-width: 64px;
  --jp-cell-prompt-font-family: 'Source Code Pro', monospace;
  --jp-cell-prompt-letter-spacing: 0px;
  --jp-cell-prompt-opacity: 1;
  --jp-cell-prompt-not-active-opacity: 0.5;
  --jp-cell-prompt-not-active-font-color: var(--md-grey-700);
  /* A custom blend of MD grey and blue 600
   * See https://meyerweb.com/eric/tools/color-blend/#546E7A:1E88E5:5:hex */
  --jp-cell-inprompt-font-color: #307fc1;
  /* A custom blend of MD grey and orange 600
   * https://meyerweb.com/eric/tools/color-blend/#546E7A:F4511E:5:hex */
  --jp-cell-outprompt-font-color: #bf5b3d;

  /* Notebook specific styles */

  --jp-notebook-padding: 10px;
  --jp-notebook-select-background: var(--jp-layout-color1);
  --jp-notebook-multiselected-color: var(--md-blue-50);

  /* The scroll padding is calculated to fill enough space at the bottom of the
  notebook to show one single-line cell (with appropriate padding) at the top
  when the notebook is scrolled all the way to the bottom. We also subtract one
  pixel so that no scrollbar appears if we have just one single-line cell in the
  notebook. This padding is to enable a 'scroll past end' feature in a notebook.
  */
  --jp-notebook-scroll-padding: calc(
    100% - var(--jp-code-font-size) * var(--jp-code-line-height) -
      var(--jp-code-padding) - var(--jp-cell-padding) - 1px
  );

  /* Rendermime styles */

  --jp-rendermime-error-background: #fdd;
  --jp-rendermime-table-row-background: var(--md-grey-100);
  --jp-rendermime-table-row-hover-background: var(--md-light-blue-50);

  /* Dialog specific styles */

  --jp-dialog-background: rgba(0, 0, 0, 0.25);

  /* Console specific styles */

  --jp-console-padding: 10px;

  /* Toolbar specific styles */

  --jp-toolbar-border-color: var(--jp-border-color1);
  --jp-toolbar-micro-height: 8px;
  --jp-toolbar-background: var(--jp-layout-color1);
  --jp-toolbar-box-shadow: 0px 0px 2px 0px rgba(0, 0, 0, 0.24);
  --jp-toolbar-header-margin: 4px 4px 0px 4px;
  --jp-toolbar-active-background: var(--md-grey-300);

  /* Input field styles */

  --jp-input-box-shadow: inset 0 0 2px var(--md-blue-300);
  --jp-input-active-background: var(--jp-layout-color1);
  --jp-input-hover-background: var(--jp-layout-color1);
  --jp-input-background: var(--md-grey-100);
  --jp-input-border-color: var(--jp-border-color1);
  --jp-input-active-border-color: var(--jp-brand-color1);
  --jp-input-active-box-shadow-color: rgba(19, 124, 189, 0.3);

  /* General editor styles */

  --jp-editor-selected-background: #d9d9d9;
  --jp-editor-selected-focused-background: #d7d4f0;
  --jp-editor-cursor-color: var(--jp-ui-font-color0);

  /* Code mirror specific styles */

  --jp-mirror-editor-keyword-color: #008000;
  --jp-mirror-editor-atom-color: #88f;
  --jp-mirror-editor-number-color: #080;
  --jp-mirror-editor-def-color: #00f;
  --jp-mirror-editor-variable-color: var(--md-grey-900);
  --jp-mirror-editor-variable-2-color: #05a;
  --jp-mirror-editor-variable-3-color: #085;
  --jp-mirror-editor-punctuation-color: #05a;
  --jp-mirror-editor-property-color: #05a;
  --jp-mirror-editor-operator-color: #aa22ff;
  --jp-mirror-editor-comment-color: #408080;
  --jp-mirror-editor-string-color: #ba2121;
  --jp-mirror-editor-string-2-color: #708;
  --jp-mirror-editor-meta-color: #aa22ff;
  --jp-mirror-editor-qualifier-color: #555;
  --jp-mirror-editor-builtin-color: #008000;
  --jp-mirror-editor-bracket-color: #997;
  --jp-mirror-editor-tag-color: #170;
  --jp-mirror-editor-attribute-color: #00c;
  --jp-mirror-editor-header-color: blue;
  --jp-mirror-editor-quote-color: #090;
  --jp-mirror-editor-link-color: #00c;
  --jp-mirror-editor-error-color: #f00;
  --jp-mirror-editor-hr-color: #999;

  /* Vega extension styles */

  --jp-vega-background: white;

  /* Sidebar-related styles */

  --jp-sidebar-min-width: 180px;

  /* Search-related styles */

  --jp-search-toggle-off-opacity: 0.5;
  --jp-search-toggle-hover-opacity: 0.8;
  --jp-search-toggle-on-opacity: 1;
  --jp-search-selected-match-background-color: rgb(245, 200, 0);
  --jp-search-selected-match-color: black;
  --jp-search-unselected-match-background-color: var(
    --jp-inverse-layout-color0
  );
  --jp-search-unselected-match-color: var(--jp-ui-inverse-font-color0);

  /* Icon colors that work well with light or dark backgrounds */
  --jp-icon-contrast-color0: var(--md-purple-600);
  --jp-icon-contrast-color1: var(--md-green-600);
  --jp-icon-contrast-color2: var(--md-pink-600);
  --jp-icon-contrast-color3: var(--md-blue-600);
}
</style>

<style type="text/css">
a.anchor-link {
   display: none;
}
.highlight  {
    margin: 0.4em;
}

/* Input area styling */
.jp-InputArea {
    overflow: hidden;
}

.jp-InputArea-editor {
    overflow: hidden;
}

@media print {
  body {
    margin: 0;
  }
}
</style>



<!-- Load mathjax -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML-full,Safe"> </script>
    <!-- MathJax configuration -->
    <script type="text/x-mathjax-config">
    init_mathjax = function() {
        if (window.MathJax) {
        // MathJax loaded
            MathJax.Hub.Config({
                TeX: {
                    equationNumbers: {
                    autoNumber: "AMS",
                    useLabelIds: true
                    }
                },
                tex2jax: {
                    inlineMath: [ ['$','$'], ["\\(","\\)"] ],
                    displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
                    processEscapes: true,
                    processEnvironments: true
                },
                displayAlign: 'center',
                CommonHTML: {
                    linebreaks: { 
                    automatic: true 
                    }
                },
                "HTML-CSS": {
                    linebreaks: { 
                    automatic: true 
                    }
                }
            });
        
            MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
        }
    }
    init_mathjax();
    </script>
    <!-- End of mathjax configuration --></head>
<body class="jp-Notebook" data-jp-theme-light="true" data-jp-theme-name="JupyterLab Light">
<div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[2]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="o">%</span><span class="k">matplotlib</span> inline
<span class="kn">import</span> <span class="nn">sys</span>
<span class="c1">#sys.path.insert(1,r&#39;G:\NewFits\Data\refnxtoolbox-master&#39;)</span>
<span class="c1">#sys.path.append(&#39;./../&#39;)</span>
<span class="c1">#import plottools, objective_processing</span>
<span class="c1">#from plottools import graph_plot, hist_plot</span>
<span class="c1">#from FreeformVFP import FreeformVFP</span>
<span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>


<span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>
<span class="kn">import</span> <span class="nn">os.path</span>
<span class="kn">import</span> <span class="nn">refnx</span><span class="o">,</span> <span class="nn">scipy</span>
<span class="c1"># the analysis module contains the curvefitting engine</span>
<span class="kn">from</span> <span class="nn">refnx.analysis</span> <span class="kn">import</span> <span class="n">CurveFitter</span><span class="p">,</span> <span class="n">Objective</span><span class="p">,</span> <span class="n">Parameter</span><span class="p">,</span> <span class="n">GlobalObjective</span><span class="p">,</span> <span class="n">process_chain</span><span class="p">,</span> <span class="n">Transform</span>
<span class="c1"># the reflect module contains functionality relevant to reflectometry</span>
<span class="kn">from</span> <span class="nn">refnx.reflect</span> <span class="kn">import</span> <span class="n">SLD</span><span class="p">,</span> <span class="n">ReflectModel</span><span class="p">,</span> <span class="n">Structure</span><span class="p">,</span> <span class="n">MixedReflectModel</span>
<span class="c1"># the ReflectDataset object will contain the data</span>
<span class="kn">from</span> <span class="nn">refnx.dataset</span> <span class="kn">import</span> <span class="n">ReflectDataset</span>
<span class="c1"># version numbers used in this analysis</span>
<span class="n">refnx</span><span class="o">.</span><span class="n">version</span><span class="o">.</span><span class="n">version</span><span class="p">,</span> <span class="n">scipy</span><span class="o">.</span><span class="n">version</span><span class="o">.</span><span class="n">version</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt">Out[2]:</div>




<div class="jp-RenderedText jp-OutputArea-output jp-OutputArea-executeResult" data-mime-type="text/plain">
<pre>(&#39;0.1.18&#39;, &#39;1.5.2&#39;)</pre>
</div>

</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[5]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">pth</span> <span class="o">=</span> <span class="s1">&#39;/mnt/StorageDevice/refnx/FinalDataD172019&#39;</span>
<span class="c1">#pth = &#39;G:\\NewFits\\Data&#39;</span>

<span class="n">data_d2oblock</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471548_471550.mft&#39;</span><span class="p">))</span>
<span class="n">data_d2oblock</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;d2osiblock&quot;</span>

<span class="n">data_d2omucins</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471667_471668.mft&#39;</span><span class="p">))</span>
<span class="n">data_d2omucins</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;d2oSilanesmucins&quot;</span>

<span class="n">data_h2omucins</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471669_471670.mft&#39;</span><span class="p">))</span>
<span class="n">data_h2omucins</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;h2oSilanesmucins&quot;</span>

<span class="n">data_d2omucinsconfined</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471678_471679_471680.mft&#39;</span><span class="p">))</span>
<span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;d2oSilanesmucinsconfined&quot;</span>

<span class="n">data_d2omucinsconfined2</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471690_471691_471692.mft&#39;</span><span class="p">))</span>
<span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;d2oSilanesmucinsconfined2&quot;</span>

<span class="n">data_d2omucinsconfined3</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471696_471697_471698.mft&#39;</span><span class="p">))</span>
<span class="n">data_d2omucinsconfined3</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;d2oSilanesmucinsconfined3&quot;</span>

<span class="n">data_hPS</span> <span class="o">=</span> <span class="n">ReflectDataset</span><span class="p">(</span><span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">pth</span><span class="p">,</span><span class="s1">&#39;471555_471556.mft&#39;</span><span class="p">))</span>
<span class="n">data_hPS</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s2">&quot;hPS&quot;</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[6]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">sldmucins</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">5.55583</span><span class="p">,</span> <span class="s1">&#39;sld mucins&#39;</span><span class="p">)</span>
<span class="n">sldmucins</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span> <span class="o">=</span> <span class="p">(</span><span class="mf">5.2</span><span class="p">,</span> <span class="mf">6.2</span><span class="p">))</span>

<span class="n">sldmucinsd2o</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">5.82485</span><span class="p">,</span> <span class="s1">&#39;sld mucins d2o&#39;</span><span class="p">)</span>
<span class="n">sldmucinsd2o</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span> <span class="o">=</span> <span class="p">(</span><span class="mf">3.7</span><span class="p">,</span> <span class="mf">5.9</span><span class="p">))</span>

<span class="n">sldmucinsh2o</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">1.8584</span><span class="p">,</span> <span class="s1">&#39;sld mucins h2o&#39;</span><span class="p">)</span>
<span class="n">sldmucinsh2o</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span> <span class="o">=</span> <span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span> <span class="o">=</span> <span class="p">(</span><span class="mf">1.0</span><span class="p">,</span> <span class="mf">6.36</span><span class="p">))</span>

<span class="n">sldhps</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">1.412</span><span class="p">,</span> <span class="s1">&#39;sld hps&#39;</span><span class="p">)</span>
<span class="n">sldhps</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span> <span class="o">=</span> <span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span> <span class="o">=</span> <span class="p">(</span><span class="mf">1.2</span><span class="p">,</span> <span class="mf">2.8</span><span class="p">))</span>

<span class="n">si</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">2.07</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">sio2</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">3.47</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">silanes</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">0.7</span><span class="o">+</span><span class="mi">0</span><span class="n">J</span><span class="p">)</span>
<span class="n">mucinsd2o</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucinsd2o</span><span class="p">)</span>
<span class="n">mucinsh2o</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucinsh2o</span><span class="p">)</span>

<span class="n">air</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">0.0</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">hps</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldhps</span><span class="p">)</span>
<span class="n">hps1</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldhps</span><span class="p">)</span>
<span class="n">hps2</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldhps</span><span class="p">)</span>
<span class="n">hps3</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldhps</span><span class="p">)</span>

<span class="n">mucinsnopocket</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucins</span><span class="p">)</span>
<span class="n">mucinspocket</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucins</span><span class="p">)</span>
<span class="n">mucinsnopocket2</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucins</span><span class="p">)</span>
<span class="n">mucinspocket2</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucins</span><span class="p">)</span>
<span class="n">mucinsnopocket3</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucins</span><span class="p">)</span>
<span class="n">mucinspocket3</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="n">sldmucins</span><span class="p">)</span>

<span class="n">pockets</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">4.94199</span> <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>
<span class="n">nopockets</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">2.61373</span> <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>

<span class="n">pockets2</span><span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">5.54975</span> <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>
<span class="n">nopockets2</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">2.60686</span> <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>

<span class="n">pockets3</span><span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">5.54975</span> <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>
<span class="n">nopockets3</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">2.60686</span>  <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>



<span class="n">Melinex</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">2.53323</span> <span class="o">+</span> <span class="mi">0</span><span class="n">J</span><span class="p">)</span>

<span class="n">Melinex</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">2.4</span><span class="p">,</span> <span class="mf">2.8</span><span class="p">))</span>
<span class="n">Melinex</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;silanes SLD&#39;</span>

<span class="n">d2o</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">6.36</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">h2o</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="o">-</span><span class="mf">0.56</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>

<span class="n">silanes</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.6</span><span class="p">,</span> <span class="mf">0.8</span><span class="p">))</span>
<span class="n">silanes</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;silanes SLD&#39;</span>

<span class="n">pockets</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4.0</span><span class="p">,</span> <span class="mf">5.0</span><span class="p">))</span>
<span class="n">pockets</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;pocket SLD 1&#39;</span>
<span class="n">nopockets</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">2.2</span><span class="p">,</span> <span class="mf">3.0</span><span class="p">))</span>
<span class="n">nopockets</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;no pocket SLD 1&#39;</span>
<span class="n">pockets2</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4.6</span><span class="p">,</span> <span class="mf">6.36</span><span class="p">))</span>
<span class="n">pockets2</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;pocket SLD 2&#39;</span>
<span class="n">nopockets2</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">2.2</span><span class="p">,</span> <span class="mf">3.0</span><span class="p">))</span>
<span class="n">nopockets2</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;no pocket SLD 2&#39;</span>
<span class="n">pockets3</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4.6</span><span class="p">,</span> <span class="mf">6.36</span><span class="p">))</span>
<span class="n">pockets3</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;pocket SLD 3&#39;</span>
<span class="n">nopockets3</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">2.2</span><span class="p">,</span> <span class="mf">3.0</span><span class="p">))</span>
<span class="n">nopockets3</span><span class="o">.</span><span class="n">real</span><span class="o">.</span><span class="n">name</span><span class="o">=</span><span class="s1">&#39;no pocket SLD 3&#39;</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[520]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#Si block D2O</span>

<span class="n">sio2_slab</span> <span class="o">=</span> <span class="n">sio2</span><span class="p">(</span><span class="mf">13.8936</span><span class="p">,</span> <span class="mf">2.19514</span><span class="p">)</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">13.8</span><span class="p">,</span> <span class="mi">14</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;sio2 thickness&#39;</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">1.9</span><span class="p">,</span> <span class="mi">8</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;sio2 roughness&#39;</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.0</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;sio2 solvation&#39;</span>

<span class="n">silanes_slab</span> <span class="o">=</span> <span class="n">silanes</span><span class="p">(</span><span class="mf">23.4608</span> <span class="p">,</span><span class="mf">10.8417</span><span class="p">)</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">22</span><span class="p">,</span> <span class="mi">24</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;silanes thickness&#39;</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">12</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;silanes/sio2 roughness&#39;</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.188672</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.01</span><span class="p">,</span> <span class="mf">0.3</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;silanes solvation&#39;</span>

<span class="n">solv_roughness1</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">1.762616</span><span class="p">,</span> <span class="s1">&#39;solvent roughness d2o&#39;</span><span class="p">)</span>
<span class="n">solv_roughness1</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">10</span><span class="p">))</span>

<span class="n">s_d2oblock</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span> <span class="n">d2o</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness1</span><span class="p">)</span>
<span class="n">s_d2oblock</span><span class="o">.</span><span class="n">contract</span> <span class="o">=</span> <span class="mf">1.5</span>

<span class="n">model_d2oblock</span> <span class="o">=</span> <span class="n">ReflectModel</span><span class="p">(</span><span class="n">s_d2oblock</span><span class="p">,</span><span class="n">scale</span><span class="o">=</span><span class="mf">0.979909</span><span class="p">,</span><span class="n">dq_type</span> <span class="o">=</span> <span class="s1">&#39;pointwise&#39;</span><span class="p">)</span>
<span class="n">model_d2oblock</span><span class="o">.</span><span class="n">scale</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.9</span><span class="p">,</span> <span class="mf">1.4</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="n">model_d2oblock</span><span class="o">.</span><span class="n">bkg</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">8.17618e-07</span><span class="p">,</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4e-8</span><span class="p">,</span> <span class="mf">6e-6</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="c1">#model_d2omucins.dq.setp(4.6689,bounds=(1, 12), vary=True)</span>

<span class="n">objective_d2oblock</span> <span class="o">=</span> <span class="n">Objective</span><span class="p">(</span><span class="n">model_d2oblock</span><span class="p">,</span> <span class="n">data_d2oblock</span><span class="p">,</span><span class="n">use_weights</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>

<span class="c1">#Si block with mucins and D2O</span>

<span class="n">mucinsd2o_slab</span> <span class="o">=</span> <span class="n">mucinsd2o</span><span class="p">(</span><span class="mf">441.043</span><span class="p">,</span> <span class="mf">6.17033</span><span class="p">)</span>
<span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">400</span><span class="p">,</span> <span class="mi">500</span><span class="p">))</span>
<span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins thickness d2o&#39;</span>
<span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">25</span><span class="p">))</span>
<span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;mucins/silanes roughness d2o&#39;</span>
<span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.924083</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.6</span><span class="p">,</span> <span class="mf">0.95</span><span class="p">))</span>
<span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins solvation d2o&#39;</span>

<span class="n">solv_roughness2</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">194.54</span><span class="p">,</span> <span class="s1">&#39;mucins/solvent roughness d2o&#39;</span><span class="p">)</span>
<span class="n">solv_roughness2</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">100</span><span class="p">,</span> <span class="mi">300</span><span class="p">))</span>

<span class="n">s_d2omucins</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span> <span class="n">mucinsd2o_slab</span> <span class="o">|</span> <span class="n">d2o</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness2</span><span class="p">)</span>
<span class="n">s_d2omucins</span><span class="o">.</span><span class="n">contract</span> <span class="o">=</span> <span class="mf">1.5</span>

<span class="n">model_d2omucins</span> <span class="o">=</span> <span class="n">ReflectModel</span><span class="p">(</span><span class="n">s_d2omucins</span><span class="p">,</span><span class="n">scale</span><span class="o">=</span><span class="mf">1.0655</span><span class="p">,</span><span class="n">dq_type</span> <span class="o">=</span> <span class="s1">&#39;pointwise&#39;</span><span class="p">)</span>
<span class="n">model_d2omucins</span><span class="o">.</span><span class="n">scale</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.9</span><span class="p">,</span> <span class="mf">1.4</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="n">model_d2omucins</span><span class="o">.</span><span class="n">bkg</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">2.32067e-06</span><span class="p">,</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4e-12</span><span class="p">,</span> <span class="mf">1e-5</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="c1">#model_d2omucins.dq.setp(4.6689,bounds=(1, 12), vary=False)</span>

<span class="n">objective_d2omucins</span> <span class="o">=</span> <span class="n">Objective</span><span class="p">(</span><span class="n">model_d2omucins</span><span class="p">,</span> <span class="n">data_d2omucins</span><span class="p">,</span><span class="n">use_weights</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>



<span class="c1">#Si block with mucins and H2O</span>

<span class="n">mucinsh2o_slab</span> <span class="o">=</span> <span class="n">mucinsh2o</span><span class="p">(</span><span class="mf">441.043</span><span class="p">,</span> <span class="mf">4.07343</span><span class="p">)</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">10</span><span class="p">,</span> <span class="mi">900</span><span class="p">))</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">constraint</span> <span class="o">=</span> <span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">thick</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins thickness d2o&#39;</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">constraint</span> <span class="o">=</span> <span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">rough</span>
<span class="c1">#mucinsh2o_slab.rough.name = name=&#39;mucins/silanes roughness&#39;</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.924083</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.01</span><span class="p">,</span> <span class="mf">1.0</span><span class="p">))</span>
<span class="n">mucinsh2o_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">constraint</span> <span class="o">=</span> <span class="n">mucinsd2o_slab</span><span class="o">.</span><span class="n">vfsolv</span>
<span class="c1">#mucinsh2o_slab.vfsolv.name = &#39;mucins h2o solvation&#39;</span>


<span class="c1">#solv_roughness2 = Parameter(97.4674 , &#39;mucins/solvent roughness&#39;)</span>
<span class="c1">#solv_roughness2.setp(vary=True, bounds=(0.01, 200))</span>

<span class="n">s_h2omucins</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span> <span class="n">mucinsh2o_slab</span> <span class="o">|</span><span class="n">h2o</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness2</span><span class="p">)</span>

<span class="n">model_h2omucins</span> <span class="o">=</span> <span class="n">ReflectModel</span><span class="p">(</span><span class="n">s_h2omucins</span><span class="p">,</span><span class="n">scale</span><span class="o">=</span><span class="mf">1.0655</span><span class="p">,</span> <span class="n">bkg</span><span class="o">=</span><span class="mf">2.03462e-06</span><span class="p">,</span> <span class="n">dq_type</span> <span class="o">=</span><span class="s1">&#39;pointwise&#39;</span><span class="p">)</span>
<span class="n">model_h2omucins</span><span class="o">.</span><span class="n">scale</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.9</span><span class="p">,</span> <span class="mf">1.4</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="n">model_h2omucins</span><span class="o">.</span><span class="n">bkg</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4e-9</span><span class="p">,</span> <span class="mf">6e-6</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="c1">#model_h2omucins.dq.setp(bounds=(4, 6), vary=True)</span>

<span class="n">objective_h2omucins</span> <span class="o">=</span> <span class="n">Objective</span><span class="p">(</span><span class="n">model_h2omucins</span><span class="p">,</span> <span class="n">data_h2omucins</span><span class="p">,</span><span class="n">use_weights</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>

<span class="n">global_objective</span> <span class="o">=</span> <span class="n">GlobalObjective</span><span class="p">([</span><span class="n">objective_d2oblock</span><span class="p">,</span><span class="n">objective_d2omucins</span><span class="p">,</span><span class="n">objective_h2omucins</span><span class="p">])</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[521]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">fitter</span> <span class="o">=</span> <span class="n">CurveFitter</span><span class="p">(</span><span class="n">global_objective</span><span class="p">)</span>
<span class="c1">#fitter.fit(&#39;least_squares&#39;);</span>
<span class="n">fitter</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="s1">&#39;differential_evolution&#39;</span><span class="p">);</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="application/vnd.jupyter.stderr">
<pre>8it [00:03,  2.67it/s]
</pre>
</div>
</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[522]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#rep = objective_processing.report()</span>
<span class="c1">#rep.process_objective(objective_d2oblock)</span>
<span class="c1">#fig, ax = objective_processing.plot_reports(rep, refl_mode=&#39;rq4&#39;)</span>
<span class="c1">#ax[2].set_xscale(&#39;log&#39;)</span>
<span class="c1">#ax[0].get_legend().remove()</span>

<span class="c1">#print(objective_d2oblock)</span>

<span class="c1"># the data</span>
<span class="n">plt</span><span class="o">.</span><span class="n">rcParams</span><span class="p">[</span><span class="s1">&#39;figure.figsize&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="p">[</span><span class="mi">15</span><span class="p">,</span> <span class="mi">12</span><span class="p">]</span>

<span class="n">plt</span><span class="o">.</span><span class="n">errorbar</span><span class="p">(</span><span class="n">data_d2oblock</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2oblock</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2oblock</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2oblock</span><span class="o">.</span><span class="n">x_err</span><span class="p">,</span>
<span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{Si\ D_2O}$&#39;</span><span class="p">,</span> <span class="n">ms</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">marker</span><span class="o">=</span><span class="s1">&#39;o&#39;</span><span class="p">,</span> <span class="n">lw</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">elinewidth</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;r&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">errorbar</span><span class="p">(</span><span class="n">data_d2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2omucins</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2omucins</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2omucins</span><span class="o">.</span><span class="n">x_err</span><span class="p">,</span>
<span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{Si\ D_2O\ mucins}$&#39;</span><span class="p">,</span> <span class="n">ms</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">marker</span><span class="o">=</span><span class="s1">&#39;o&#39;</span><span class="p">,</span> <span class="n">lw</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">elinewidth</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;b&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">errorbar</span><span class="p">(</span><span class="n">data_h2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_h2omucins</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_h2omucins</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_h2omucins</span><span class="o">.</span><span class="n">x_err</span><span class="p">,</span>
<span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{Si\ H_2O\ mucins}$&#39;</span><span class="p">,</span> <span class="n">ms</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">marker</span><span class="o">=</span><span class="s1">&#39;o&#39;</span><span class="p">,</span> <span class="n">lw</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">elinewidth</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;g&#39;</span><span class="p">)</span>
<span class="c1"># the median of the posterior</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">data_d2oblock</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">objective_d2oblock</span><span class="o">.</span><span class="n">generative</span><span class="p">(),</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;r&#39;</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">20</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">data_d2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">objective_d2omucins</span><span class="o">.</span><span class="n">generative</span><span class="p">(),</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;b&#39;</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">20</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">data_h2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">objective_h2omucins</span><span class="o">.</span><span class="n">generative</span><span class="p">(),</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;g&#39;</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">20</span><span class="p">)</span>
<span class="c1"># plot the spread of the fits for the different datasets</span>
<span class="c1">#gen = objective_d2oblock.pgen(500)</span>
<span class="c1">#save_parsfucins = np.copy(objective_d2omucins.parameters)</span>
<span class="c1">#for i in range(100):</span>
<span class="c1">#    objective_d2omucins.setp(next(gen))</span>
<span class="c1">#    plt.plot(data_d2omucins.x, objective_d2omucins.generative(),</span>
<span class="c1">#    color=&#39;k&#39;, alpha=0.02, zorder=10)</span>
<span class="c1"># put back the saved parameters</span>
<span class="c1">#objective_d2omucins.setp(save_parsfucins)</span>
<span class="n">ax</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">gca</span><span class="p">()</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xscale</span><span class="p">(</span><span class="s2">&quot;log&quot;</span><span class="p">,</span> <span class="n">nonpositive</span><span class="o">=</span><span class="s1">&#39;clip&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_yscale</span><span class="p">(</span><span class="s2">&quot;log&quot;</span><span class="p">,</span> <span class="n">nonpositive</span><span class="o">=</span><span class="s1">&#39;clip&#39;</span><span class="p">)</span>
<span class="c1"># ax.text(-0.04, 1e-11, &#39;a)&#39;)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">legend</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">yscale</span><span class="p">(</span><span class="s1">&#39;log&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylabel</span><span class="p">(</span><span class="s1">&#39;Reflectivity&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlabel</span><span class="p">(</span><span class="s1">&#39;Q /$\AA^{-1}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylim</span><span class="p">(</span><span class="mf">1e-8</span><span class="p">,</span> <span class="mi">2</span><span class="p">);</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlim</span><span class="p">(</span><span class="mf">0.004</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">);</span>
<span class="n">plt</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;Freemucinsd2o.pdf&#39;</span><span class="p">)</span>

<span class="nb">print</span><span class="p">(</span><span class="n">global_objective</span><span class="p">)</span>
<span class="c1">#np.savetxt(&#39;fitmucinsd2o.txt&#39;,save_parsfucins)</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="text/plain">
<pre>_______________________________________________________________________________

--Global Objective--
________________________________________________________________________________
Objective - 2289404052872
Dataset = d2osiblock
datapoints = 312
chi2 = 353.3083170711377
Weighted = True
Transform = None
________________________________________________________________________________
Parameters:       &#39;&#39;       
________________________________________________________________________________
Parameters: &#39;instrument parameters&#39;
&lt;Parameter:    &#39;scale&#39;    , value=0.979909 (fixed)  , bounds=[0.9, 1.4]&gt;
&lt;Parameter:     &#39;bkg&#39;     , value=8.17618e-07 (fixed)  , bounds=[4e-08, 6e-06]&gt;
&lt;Parameter:&#39;dq - resolution&#39;, value=5 (fixed)  , bounds=[-inf, inf]&gt;
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=6.36 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;solvent roughness d2o&#39;, value=1 +/- 0.364, bounds=[1.0, 10.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;


________________________________________________________________________________
Objective - 2289403889032
Dataset = d2oSilanesmucins
datapoints = 309
chi2 = 931.8099383462879
Weighted = True
Transform = None
________________________________________________________________________________
Parameters:       &#39;&#39;       
________________________________________________________________________________
Parameters: &#39;instrument parameters&#39;
&lt;Parameter:    &#39;scale&#39;    , value=1.0655 (fixed)  , bounds=[0.9, 1.4]&gt;
&lt;Parameter:     &#39;bkg&#39;     , value=2.32067e-06 (fixed)  , bounds=[4e-12, 1e-05]&gt;
&lt;Parameter:&#39;dq - resolution&#39;, value=5 (fixed)  , bounds=[-inf, inf]&gt;
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;mucins thickness d2o&#39;, value=441.043 (fixed)  , bounds=[400.0, 500.0]&gt;
&lt;Parameter:&#39;sld mucins d2o&#39;, value=5.55583 +/- 0.0395, bounds=[3.7, 5.9]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;mucins/silanes roughness d2o&#39;, value=6.17033 (fixed)  , bounds=[0.001, 25.0]&gt;
&lt;Parameter:&#39;mucins solvation d2o&#39;, value=0.924083 (fixed)  , bounds=[0.6, 0.95]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=6.36 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;mucins/solvent roughness d2o&#39;, value=194.54 (fixed)  , bounds=[100.0, 300.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;


________________________________________________________________________________
Objective - 2289427497992
Dataset = h2oSilanesmucins
datapoints = 309
chi2 = 422.04634749174767
Weighted = True
Transform = None
________________________________________________________________________________
Parameters:       &#39;&#39;       
________________________________________________________________________________
Parameters: &#39;instrument parameters&#39;
&lt;Parameter:    &#39;scale&#39;    , value=1.0655 (fixed)  , bounds=[0.9, 1.4]&gt;
&lt;Parameter:     &#39;bkg&#39;     , value=2.03462e-06 (fixed)  , bounds=[4e-09, 6e-06]&gt;
&lt;Parameter:&#39;dq - resolution&#39;, value=5 (fixed)  , bounds=[-inf, inf]&gt;
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;mucins thickness d2o&#39;, value=441.043          , bounds=[10.0, 900.0], constraint=&lt;Parameter:&#39;mucins thickness d2o&#39;, value=441.043 (fixed)  , bounds=[400.0, 500.0]&gt;&gt;
&lt;Parameter:&#39;sld mucins h2o&#39;, value=1.8584 +/- 0.0985, bounds=[1.0, 6.36]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=6.17033          , bounds=[0.001, 30.0], constraint=&lt;Parameter:&#39;mucins/silanes roughness d2o&#39;, value=6.17033 (fixed)  , bounds=[0.001, 25.0]&gt;&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0.924083          , bounds=[0.01, 1.0], constraint=&lt;Parameter:&#39;mucins solvation d2o&#39;, value=0.924083 (fixed)  , bounds=[0.6, 0.95]&gt;&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=-0.56 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;mucins/solvent roughness d2o&#39;, value=194.54 (fixed)  , bounds=[100.0, 300.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;


</pre>
</div>
</div>

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>




<div class="jp-RenderedImage jp-OutputArea-output ">
<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA4EAAALCCAYAAABgEda6AAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjMuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8vihELAAAACXBIWXMAAAsTAAALEwEAmpwYAADCCElEQVR4nOzde3zbVf3H8df5ps3u96bAuIylHcuFweQmoj835apSQQRxCIgKigOGPxRvKE1RvPxAdAOG4g1QHBcRNPMCzLkpyEXukLSjbbaxyaXdVi67tk3O749vkyZtel3vfT99xCbf60nasX72OefzMdZaREREREREZHRwBnsAIiIiIiIiMnAUBIqIiIiIiIwiCgJFRERERERGEQWBIiIiIiIio4iCQBERERERkVFEQaCIiIiIiMgoUjDYA+gPRUVF9uCDDx7sYYiIiIiIiAyKp59+eou11pdv34gMAg8++GCeeuqpwR6GiIiIiIjIoDDGbOxon6aDioiIiIiIjCIKAkVEREREREaR0RUERiKDPQIREREREZFBZay1gz2GPnfUUUfZvGsCjYER+H5FRERERPpCU1MTmzdvZvfu3YM9FOmmsWPHcsABB1BYWJiz3RjztLX2qHznjMjCMCIiIiIi0nObN29m0qRJHHzwwRhjBns40gVrLVu3bmXz5s3Mnj272+eNrumgncg3U1SzR0VERERkNNm9ezczZsxQADhMGGOYMWNGjzO3CgJbVFR0b5uIiIiIyEimAHB46c33a3QEgYkEhMPu83DYfd3bw5UeFBERERGRYWx0BIFlZVBV5T6vqnJft8gX8HVyuNKDIiIiIiJtKVEyrIyOILCyElIp93kqBZWVmZ/TfAFfnsNFRERERKQjfZwoufbaawmHwxx22GHMnz+fJ554AoDjjjsu7/Eej4f58+cTDoc5/PDDueGGG0ilf6HPY/PmzZx22mnMmTOHkpISLr/8chobG/v0PQxloyIIrCsKgtPyVh0HgsHMz2m7gC9uCbY/vMdTStP0jyIiIiIiIt332GOPsXLlSp555hleeOEFVq1axYEHHgjAv//977znjBs3jueee45YLMbDDz/MX/7yFyo6CEyttZxxxhmcfvrpVFdX8/LLL7N9+3auuuqqfntPQ82oCAKPrY9CIOC+CAQgGs3saxfwEeeEE9ofXndsZ3NEO5b+2RvOwWA6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gz0yERERkRFp5AeBZWU5LR66vRZQmUGRgdFSInQcu9nMgXhIum1cVl0+2CMTERERGZFGZBD46sSsF5WVOS0e5rKuexdRZlBkwCw9IcpbTOZV9ucwnieFR/0CRUREhpHImshgD0F6YEQGga9NynoRDGZaPCRMCQHcrGDY90bPak8oMyjSby5f6ic58yAA5vIyjoP6BYqIiAwjFWsr+vR61157LeFwmMMOO4z58+fzxBNPAHDcccflPX7ixIk5r2+77TYuvfTSDq+/efNmTjvtNObMmUNJSQmXX345jY2NffcGuqGj9zIQRmQQmCOrxUNZ4d+opQSAqq3FPas9od9GRfrVjXdMBWAGW9WNRUREZBR77LHHWLlyJc888wwvvPACq1at4sADDwTg3//+915f31rLGWecwemnn051dTUvv/wy27dv56qrrtrra/dEX7yX3hr5QWBWi4fKplLAAO4MUU03Exk6FiyAcZ492LHjiL2YUjcWERGRYSDRkCC83C3CGF4eJtGw922eXnvtNYqKihgzZgwARUVFzJw5E2if8euN1atXM3bsWD7zmc8A4PF4+PGPf8yvfvUrdu7cmXPshg0bCAQCXHjhhRx66KF86lOfYtWqVbz3ve9lzpw5PPnkk5njDj300Mx5119/PZGWQOOOO+7gsMMO4/DDD+e8887LHJN+Lxs2bCAYDHLRRRcRDoc56aST2LVrFzt27OAjH/kIhx9+OIceeih33333Xr/3tBEbBOb7IQwGW59rupnI0FJQAO+d/Rprdh8LhYXqFygiIjIMlK0oo2qLu9yqaksVZSv2vs3TSSedxKZNmzjkkENYvHgxa9eu7fKcXbt2MX/+/Mzj6quv7vDYWCzGkUcembNt8uTJHHTQQdTU1LQ7vqamhssvv5wXXniBqqoqfve73/HII49w/fXX873vfa/TccViMa699lpWr17N888/z9Kl+aufV1dXc8kllxCLxZg6dSr33Xcff/vb35g5cybPP/88L730EqecckqXn0N3jdggMN8PYfb0Mk03Exl6PrDtPl7kMLakpqlfoIiIyBAVWRPBVBhMhSFeHydl3SKMKZsiXh/P7OttsZiJEyfy9NNPc+utt+Lz+Tj77LO57bbbOj1n3LhxPPfcc5nHNddcA8ADDzzARRddxGmnncZDDz0EuNNBjTHtrtHR9tmzZzNv3jwcxyEcDnP88cdjjGHevHls2LCh03GtXr2aM888k6KiIgCmT5+e97jZs2czf/58AI488kg2bNjAvHnzWLVqFV/72tf417/+xZQpUzq9V0+M2CAwZVNU1lfmbMueXhaLoelmIkNJJMLCbfcBsJYFmrMtIiIyREUWRrDlFltuCflCOMYNKRzjEPKFMvsiCyO9vofH42HhwoVUVFRw0003cd999/XqOqeffjo///nPue222zLTKcPhME899VTOcW+//TabNm2ipKSk3TXS01IBHMfJvHYch+bmZgAKCgpIpTsSALt37wY6Diw7u4fH46G5uZlDDjmEp59+mnnz5vGNb3wjE9j2hREbBDrGIegLqlytyHARiXB0cAfj2cHfOIUx7CLsrSZxfmSwRyYiIiIdiC6KEihyizAGigJEF+39VLt169ZRXV2def3cc88xa9asvbrmd7/7XS655BIAjj/+eHbu3Mkdd9wBQDKZ5Mtf/jIXXHAB48eP79X199lnH+rq6ti6dSt79uxh5cqVmXvdc889bN26FYBt27Z1+5qvvvoq48eP59xzz+UrX/kKzzzzTK/Gls+IDQLTP4Rty9Wq04PI0FW48n6KzBYe4z28myepai7RjFAREZEhzD/NT2yxW4QxtjiGf9reT7Xbvn07n/70pwmFQhx22GHE4/FMkZWestbyta99jQ996EMcccQRABhjuP/++7n33nuZM2cOhxxyCGPHju1yfV9nCgsLufrqq3n3u9/NqaeeSqClO0E4HOaqq65iwYIFHH744VxxxRXdvuaLL77IMcccw/z587n22mv51re+1evxtWWstX12saHCzDTWvuq+L1NhsOUWjIGW95r1VESGkEgEvl+xh0bGEOFqIlyDMe7MUBEREel/lZWVBLOrKXZT5nfuIWbZsmXcfvvtHH300cyfP5+LL754sIfUL/J934wxT1trj8p3/IjNBGoaqMjwE4nAzIPdOfHTaMAxll78PSQiIiIDrHzB0Jxut2TJEp5++ml++tOfjtgAsDdGbBDYdhpoVxQ0igwNDz4I48xunuTdfNTeT7TxZLWKEBERGeL2pgiMDLwRGwT2VE+DRhHpH4ccAosmreRPfJS7WIQ/sUqtIkRERET60IgOAsPLw5mviWlZO/QvFSJDVyTCmW//kneYzMOcqFYRIiIiIn1sRAaB+03aD3Abxqe/li3KOmBhz6eKarqoyACJRDg++BqTeYvfcybH8LhaRYiIiAygkVg4ciTrzfdrRAaB9TvqAbdhfPprpS//sYmGRG7GsKF17VE68KtYW9Gr6aIKHEV6x7vyD3hNM3/kND7BPWoVISIiMkDGjh3L1q1bFQgOE9Zatm7dytixY3t0XkE/jWdQNdvmnNeOcQjU5a8xX7aiLDdjuKIs0+tkb9cJVqyt0CJZkV6I3OFnS8vfPfvzX5xUM5WVhYM7KBERkVHggAMOYPPmzdTX1w/2UKSbxo4dywEHHNCjc0ZkEEibf7iYMW4GJyTy/yDH6+OZ5ymborK+Mme/CsaIDLxIBO6+G2qrGnmAj/Eh8yC1wVMHe1giIiIjXmFhIbNnzx7sYUg/G5HTQTHuF8e4b6/uyjqWHdvx4enjHOMQ9HXclGzOjXNypot2pLMppiLSPX/+M0yZAg9wOmfYe1gXTxKe06huESIiIiJ7aUQGgWML3DmxgaJAt45PHxcoChBdFO3wuJptNZSt6HphUr4ppiLSM34/rJ7xCRoZw1aKKOYNqmo8WhsoIiIispdGZBC4/6H7A2TW9qULtHSUoUsfF1scwz/NT6IhwZwb5+S9dtvpomnZRWAq6ytzitLE6+MqEiPSU5EI4cSfGMsufs1n+Ry/IIWHyrgWqouIiIjsDTMSK/+YCmMBjj3gWB7f/Hhme8gXompLVSZA83q8NCYbCflCxOvj2HL3swgvD+esFcwW8oUyQWObe2bOL76umK27tpKyKbcoTVEg7zki0oVwmJnxh3iN/XmAj3IW9zLHu4lYpeOmCkVEREQkL2PM09bao/LtG5GZwEP+fQgAJr04sEW8Pp4JAAEak41Aaz9BgIW3LewwACydXpozXTSd3Wub5avfWd/tKabZlC0UaSMa5W/+S/Cyh7tYxMX8lGjTKWhOqIiIiEjvjcggcNLWSQD8+3P/ztk+xjMm7/HpwDCyJsLajWsJ+UJ5j6u+rBr/tNbsQ8XaChINiUwF0c6mmLaVL+BTJVKRNvx+DjtvPktYxr2cxVf4EX5bC/G4W0JURERERHpsRAaB2ZYcsyTz/MApBzJ1zNQOj00HYe878H2UTi/t1vWzi75U1lcSvNmtLppee9jVvUDVREU6FYlw2ZwHAbiRyziaJwh7q0mcHxnccYmIiIgMUyM+CFy1flXmeaIhwZt73sy8Dha5Adu+E/elZFpJZvutz9yKb7yv3bWyC8wUX1cM5PYZtNicKab5rpEv4FM1UZHOHfS3W5li3uFWPs8X+BlVzSWaESoiIiLSSyM2CCxfUA7kVvPMXg8IEL/EDeBe+/Jr1CypydmXtMl216xYW8Hr21+nbEUZ9TvzN5/PvteWnVvabW8b8B37i2Nz1irma1gvMtpF7vCzzU7jbaawjRnsn9rEungzahooIiIi0nMjNgiMLIwA5DR/TzeFT0tn5CJrIu2mYK74+IrM88P3OZwJhRMA2O9H+3VYOCb7HgbTrvF8ZE2kXcC3ZecWQr5QtxrWq3CMjFaRCIRCMJZd3MAVXMKNzGWdCsSIiIiI9MKIDQLBnXqZnp4JtCvQks7IVaytyKzlS8uekvncxc/x9jfeBuD7x3+f8YXjO7xnuiqoxbarCpoOTNsGfNFF0W5VE1XhGBnNoicsZR9e5w32xcHyZz6sAjEiIiIivTCig8CyFWU5GT6vx5uzP1+7iLS22b504Pb1932dF7/4Yof3nDp2aub57Kmzc/alM3ltAz7/NH+mmuhZobMywWpXTe5FRhP/0svZEDqV9/BvfswV/IbzVCBGREREpBdGdBBYWV+ZE+i1XWvXdnpoNq/Hy5JjlmTWFqalC7nks+SYJby95+3M6zHfHcP595/P32r+RmV9ZbtMXtv2EeULynOOST/PXkdYWV+pwjEyekWjvFowi9eYiZcm3miaphmhIiIiIj00ooPAoC+YE+i1XWuXzshB+8byzalmVq1fRWRhJCcTF7w5mNNcPtvSDy3NCTqbUk389oXf8qE7P5TTMiJf4ZfImkhmumhb2cGsxRKvjysbKKNS5A4/G5v3ZwLv8CO+zEX2Z8TjlvCcRtWIEREREemmER0EZq+1g9Ypn+nWDbHFsUymr22AmLIp4vVxImsiOZm4xmRjuyqj2doGeBbLLz/6Syw2ZxtA6bJSNry5AXCzfulAcc6Nc5hz4xzAnf5ZMr2EtpQNlNEoXSBmJ+OppxiLYQZbqarxKCMoIiIi0k0jOgjMXmsHZLJnW3dtzWxLZ9+yjwN3qmjIFyKyMNJuWmk+IV8IyM0+pq/xo8d+1O74MZ4x1DbUMnvpbN71s3cBrQFkzbYaarbV5GxrKx2gQm7VUFUQlZEuesJSglQxibe4mUu5jGWk8KhGjIiIiEg3jeggsK3s1gxdya7SmR3YGUymwIzHeAA3s5g+Nl+lz3yBXLpH4Tfe9w3WbVkHkJMtTLNYarfV5jSeTweX5x9+PuHl4UwWMdGQUAVRGfH8Sy8nFvoExdSxk/Fs4iAOJkHIW0PkfM0JFREREenKqAoCszN00L7qZrbsoi3ZgV3QF6TyEjeoa766GYC6K+syx2ZnH9PXaDvVtHR6aWY6548e+xF7kns6HbfH8VC/sz4TCKaDy2N/cWxOwZh0m4vw8jCX//Xybn8uyh7KsBON8lDppezPf/kVn2UJy4g2naK+gSIiIiLdMKqCwOwMHeRW3eyo2AvkD+x6Il/fv87WGBaYgpzXzSk32GxKNQHw3Bee447n76B+Z31OwZj0mseqLVUse3IZ0HmAlw6CszOJIsOC34//U+/hBQ6jiC3cx5nMtrXqGygiIiLSDaMqCMwO5CC36mZ3poj2Rr6WErXbavPeL52hTJF/LG/ufhOA4uuLWfrE0g7vmb52OsDrSHYmsWpLlYrNyPASiTA1tD/f55s8yvv4Hee4VWMUBIqIiIh0alQFgW21LeLSl9JVR7OzjeCuH2zbuiK9xnDGuBlA1wGptTYTEKYVOAXtjks3vO8oy5edSUzZVLu1i5omKkNeNMqCks0cyCtcyXUsarhZrSJEREREujBqgsC2Td+hfRGXjo7ryXUzFTsXRoisiRCvj+cEdVt2biG6KJozpfSgKQcB7tpC33hfhwFpuiBMdkP6tPSU0XzaNphPNCQyLSiyr100viizX9NEZVjw+zltzN+ox8drzKTwtY2UBapRJCgiIiLSsVETBLZtxF6+oLzdWr/yBeUdNmzPXKdNdqzt8dnTL9P7siuLBn1B/NP8mewfkBNkPX7h45mANN12Iv01XRAmu9CMYxyKxhXxwNkPdDhmi83J8h37i2MzLSjSUjZF/c56wsvDnPzbkzVNVIaFSMRdBribcZRSzd18kjlNL1J3rH5mRURERDoyaoLAtvIFe10FgNAa5PUkY5hdWTRdJCY7KEtnCiNrIu0C03xf0wVgwC1a88RFT3Ba4LRMsJjP9HHTidXFiKyJUL+zvsPjqrZUUbOtptNpoiJDRSQCIV8dDklqKKGAZhqYwUfqf0VixtHKCIqIiIjkMWqDwM50FOBlb+9OwJiWr7Jo22xeujF9R/fMXmOYnTn0erw5rSzSgeD0sdMz2ycUTqBhdwOH3nIoP3jkB52Ote16RMc47Vpc5KP1gzJYoo8XE/CuBwxF1PNPFlBCLWXbblPLCBEREZE8RkUQ2NN1fh0FeD0J/LoSXRRt1/cvW3rM6Xt2tMYwO0uXnUXc+rWt1C6ppXxBOdu/uZ1Xr3iVfSfum7cnYduWFNn80/ycMPuEnG3ZAZ/WD8pg8/vh8xc7eEjyCgdTQg0PcQopDKl4paqFioiIiLQxKoLAvgzeuquzRvTgBld1V9YB+XsPdjRdNeQL5VQ0TWfpsgOzfBnLfSbuwxvb38g71mbbcVGZV956hWVPLiO8PMyaDWvaBXxtey1q/aAMhsuX+pkbKsAhyVZmsIMJ7MMbzONFEt+5E8JhTQ0VERERaTEqgsCB0Dboa1tcpa+0rWiaziDmK0jTVtH4IgwGaC0oM3Xs1E7vl92APvs9pSuOtu21GK+Pa2qoDIpoFAKlSd5kKvN5lrV8gGlsoyz1AFRVaWqoiIiISAsFgX2kbUasbXEV6Jt1c/kKx2QHn/mmY17+18sJLw9Tv7MeiwXcAPKJi57grd1vdeu+KZuiMdmYeS8WS7w+zrSx09q1tLg3fq+mhcqA8/vhrE95AcN/OJpSqqlhDpvYH1IpqFSBIxERERFQENhtbTN92UFO27V66a/Z0za7KvzSU/ma0Xc0HXPZk8vaVfhM9yrMblxvMBQ6hR3e02M8ebenM5NpmhYqgyUSgVAIjONhLLuox8ehxAnzEomSEwd7eCIiIiJDgoLAbuos2Mq3Vq90emneaZtt9Xa9Yvq8ttMxs4O97MA1nQHMfj+QO7006AtSdWnHU1eTNtlu27bd26jeWp2zLT0tVIViZDBEoxAIwEvM4wOs5jGOI4WhjPx/BkVERERGGwWB3dTV2re2a/UePPfBvK0h+lp2Jq9tO4fswLWt9PjT00vLF5Rnxlm+oBzfeF/ONM9gUZC/n//3vNdqSjVl7p9NGUEZDH4/nHUWgOEJ3s0c1rGNGdTWNFNXrAIxIiIiIqMuCMw3/bI7UzLbBlttp3e2XavXX0FfWx0VioHcwDVbvvG3ff74hY9nrusb72P5R5Zz2V8vc8/v4McmXXQmTYViZLCkp4XudCazD2+wjekcSpzirSoQIyIiIjLqgsCOWi90pbNgazB1Fny2Xe+X1p3xZ1+37so6LvnLJa1ZRZN/fWC+6aK+8b6czzc9RdWpcDRdVPpVNAqBGXX8m+P4GH/gaY7iyNTjJOK71DtQRERERrVRFwT21mBl+vZG2/V+S45ZAvRs/OULyvMWvknaJF6Pt8vz63fWYyoMc26ck9NX0GKp2lLFsb84tpfvTqRzfj/E6ooJeDewhgUcyVNUcwgfLnhIQaCIiIiMagoCR4h8U1rbBq5LP7S0x9eNLIzkLXwT8oWovKQS33gfACFfiNLppe3WBabVbKvhpN+c1G5tZf3OeoqvK1ZGUPpFJALxxhLq2ZeZ/JdmCtnT7KHJMwYKCtREXkREREYlBYEjRGdTWrMDxN62pMg3HdY/zU/dlXWZojIPnvtgu3YR2WobahlXOK7d2sGtu7aqgIz0C3dtoMFxIMpH+TS3s4HZzE89RSJ5kJrIi4iIyKikILCf9Tbo6ksdFYDpCf80P2eFzsqpItr2mtmVRkO+ULtrFE8oZtrYae3aVbRtbSHSl9JrAwHu5zRO5CGqmcsp/MVtIh+Pa3qoiIiIjCoKAvtY26Cvt0HXUFSxtqJb7yeyMEJ0UZTS6aWZbaXTS3nsc4+x8UsbWXToIrxO7npCxziULivFqXA0PVT6lN8PZy0uxkOSN5jJGHZTxBYamME7THTLiCoIFBERkVFEQWAfG0lBX1p20/nuVvT0T/NTfVk15QvKKV9QTvVlbkP5w356GCteWsFBUw+iwBRkjk/aJLUNtVgs9TvrCd4czLmP2kzI3ohEYG5pCockK/ko53M7W5lBKdWEG5/VskAREREZVYy1tuujhpmjjjrKPvXUU/1ybVNhsOUj7zPrTHh5mKotVaRsCsc4BIoCmYIz3RVZE+He+L0518nXwzDb7KmzWXX+KspWlBGvjxPyhTJrEUV6KpFwl/+tWwfTnTc5q+m3LOdSxrCHktAYYj37kRYREREZ0owxT1trj8q3T5lA6VS+9hA9Wb+XziJWrK1odx2gw2qiAOvfXM/cm+Zm7le1pUoFZKTX/H6IxeBb34L6pqnECXEMT+CQJB637DujURlBERERGRUUBEqn8rWHCPqC3T4/3RewLcc4lE4vzVQTLZ1eSqFT2O645lRzppBMyqaI18cJLw9jKgxjvjtGawilx9yKofBPFnAiD1FAM1N4i7ptHhUKFRERkVFBQaB0KV97iO7K7guYrcAp4MFzH8xUE62+rJovHvXFTN/BzsTr4wA0JhszawhLlpXgVDjdXrMoo1v0hKUEqOJHfJmrqeAtpjKOXSoUKiIiIqOCgkDpUtum8z1Zkxf0Bdv1BQRoSjZlrnP+4ecTXh5m2ZPLOryOweTNFGaz2Jwpo+mpqE6Fw5wb5zDnxjkKFAUA/9LLiYU+gZ/1/JrPsJib2clEDEnuvVf940VERGRkUxDYQ0Oh799g6c17jy6Ktps+2nZKafaU0a27tuL1eDPTTw2G0umlFI0voinV1OX90lNGnQqH4M1BqrZUYbHUbKuhZltNu0BRRrFolGjpFYBhO+M5lscopJnKeJKykhiEw4oGRUREZERSENhDI7EFRHf15r1nZxHTDeSzp5TmKzzTmGzMZAktlo1vbmTrrq09uq/F0phszDsVNR0oqu3EKOf3469+kLPKw9zBBZzDnUziHcaxmyrmkqqsQosERUREZCRSECgDonxBeWb9X/aU0nyFZ0K+EF6PNzONtCnV1C6YG184nn98+h+ZNYTZPQe7kr7HaA7opVUkAqHSJr7Nd/g/rmQX4yigmUPtCyTiu7RIUEREREYc9QmUQZdoSLTrBVi6rDRTFTRbR30Kw8vDVNZXtjvHYLBYSqaVYIyhZlsNXo+XpmQTQV9QfQcFgMSckymruYHx7OBkHuRavo0hSZAqYqFPQDTq9pgQERERGSbUJ1CGtPSU0ewsYdAXzFkX6PV4gY6rk3YWAAJsfGsjNdtqKHQKM20ntDZQAIhE8Nc8RIxDeZYjeIdJnMDDOFhNCxUREZERSUGgDBnZ0zPTbSkMhqAvSOUlle2mkmbLDhod47jTSU1rVdLmVDOQO7VUawMFaG0c6DjMZR03cSmfZAX78joT2MnbdhLE4+A4KhYjIiIiI4Kmg8qIkD2l1DfeR/3O+m6dV+AUkEwlNTV0tEskoKyMRFUjZQV/5dXG6dzKRZzDXRTQzGwSrKQMv7MRAgGIxbq+poiIiMgg0nRQGfGyp5TWXVmXU2wmW9t+g+mpoZX1lZz0m5MILw9TcE2BegmONn4/xGL4k9XE9pTSUPsmdxZ8lv/jq+xmHOuYSxlRSKVQR3kREREZ7hQEyoiSnlKank4KUDq9lNLppZmppVWXVmHLczPgFkttQy3x+jhJm9R6wVEucoefPzZ/hKc4ks/yS1IUUEWAlHHcqaMKAkVERGQYUxAoI1J2ZrD6smqqL6smVZ7KWVPYUbYQWtcLyuiUXiZ4l3Muh/Mcx/IYDile8syHdeu0NlBERESGNQWBMqJ11Asw0ZDosJl8NlNhCN4c1NTQUSgadZf/fd25nu+M+S4+6jm9+V62JqdAlSqGioiIyPClIFBGpbIVZZ0Gdk7WH42qLVWULCuh8JpCTIXBVBjm3DhHgeEI17JMkJeqx/C9CdfyCz7Hf9mf0/gjzSmjtYEiIiIybCkIlFEnsiZCvD6eNwtoMHiMhxTt9zXb5szzmm01mTWDiYYE4eVhnApHBWVGoLIyWPvmfK7lW9zEJTzK+/gCPyPsrabguxHNDBUREZFhRy0iZFQKLw9TtaWKlE25FUM9hTQlmyj0uM3ku5ommuYb72PK2CkkGhLtzimdXsqD5z6othPDWCQCFRWtry/iVsaymxtZgiGFxcEhSaA0SazaO2jjFBEREWlLLSJE2sjXjD5VnqIp2dTtABCgfmc9Ndtq8p6TnS2U4SmrjzwAP+fzzKWKhfwDg/sPaCk8rKtxNDNUREREhg0FgTIqpauHtq0YGvQFMxVDHeNQOr2UkC8EkNNfsLvi9XEiayJ9Nm4ZeNEozJjR+vp/+QlX8CMOYDMFNGFIMpd1RJYXa16oiIiIDAsKAkWyZGcIA0UBHjz3QWKLY9hyS+O3G7HlttPWEm2FfKEOK5TK8OD3Q11da0awCS8XcSu/4HOMZTcT2MldnA1bt6piqIiIiAwLCgJFsnSUIcwWXRRlxrgZec7OdeDkA4kuivbHMGUQpFtGGAMzSosITn6NezmLnYznKr5HMoUqhoqIiMiwoCBQpIf80/zUXVmXkxF0jEPIF6LuK3V8433fYKJ3Ipve3sRhtxzWYVsJVRUdXtItI1IpiFV7OeCtGKccGGcZS4jyUb7KdW66UEGgiIiIDHEKAkV6qe3U0eiiKL4JPr53/PdYe8FaHOOwo2lH5vh0oZh08FeyrIR4fRyLJV4fp2RZiYLBYSZxxyPs8EzhMpZxA1dQXPVPCkwz4TE1JNa8MtjDExEREclLLSJE+kHxdcXU76zv8XmOcQgUBYgtjgFutrBsRRmV9ZUEfUGii6JqOTGEhMNQVWm50S7mz5zKX/lQa9sI73pie0oHe4giIiIySqlFhMgAiqyJ9CoABEjZFPH6OMXXFWcCwKotVVgsVVuq1HJiCIlE3CWAKWu4nGV8gZ8SJoaHZrdtROPBRBauGexhioiIiLSjIFCkj0UWRgj5QhhMr6+xddfWTAYw3YMwHSCq5cTQkN1DsJlCPs3tLGUJM9hKAY2UUEOk/hK1jRAREZEhR0GgSD+ILooS9AUB8Hq83W4pkZYO+Eqml7QrPqOWE0NHumKoxwNFsyZybUGE33AuHlIkKCEYv5dE4MNQUODOHVVAKCIiIkOAgkCRfpBuNVG+oJzKSyoJFAXaHdPdTGHb4jMydKQrhjY3Q/UGL3+/6h/8hvO4nU+TxMPLzOXUpvsgmYSqKvURFBERkSFBQaBIP4osjGQCQltuqV1Sm5kqGvQFKZ1e2mmWsGZbDWeFzuq0b6EMHREirOAcVnIqP+Z/SeGhigAW3N4S6iMoIiIiQ4Cqg4oMorbVP3c07mDjWxsz+x3jkLIpQr6QKoMOE+E5jVTVePgaP6QRLz/iKxiSBKkiWnoF/uoHB3uIIiIiMgp0Vh1UQaDIEJJoSHDyb0+mZltNznaDYfa02TjGyewrnV7Kg+c+mBMYqqXE4Esk3FmflXHLzeYS/mn/h7tYhCFFsLSZWLV3sIcoIiIio4BaRIgME/5pfrweb7v1ghZLoiGRExymm89nU0uJwZdeJ3h1ueFSeyMf4w98gNUAxGu8FBerPoyIiIgMLgWBIkNMZX0llu5l6LN7CkbWRIjXx3NaSlTWV/bnUKUTkQgEQh4+Y+7gq/ww00NwyxbVhxEREZHBpSBQZIgJ+oI9ailRv7OekmUlXLP2mpx2FI5xMm0qZHBEozChaByf4nfcyKXswxs4tol4HIxRjRgREREZHAoCRYaY6KJopi1EyBfKqSjaGYulMdlIyqbwGI9aSgwBfj/U1cG+pZP4LL/iTs5hAjsZX9jI1q0KAkVERGRwKAgUGWLSLSWy20Kkm8+nA8PS6aWdXmNu0Vy1lBhCopQxjt18heu5h7NobrKUzXiUWMlHCc9pxHHUS15EREQGjqqDigxDiYYEwZuDNCYbOz1OrSWGgEgEKipyNv2ej/MJ7mEi77CdiVg8OA4EAm5RGREREZG9peqgIiOMf5qfyksqCflCnR4Xr49TsqyE8PIwiQalmQZFJAKhEDjuf24TzOZJjuZnfIF3mAIt03zdXvJWU0RFRESk3ykTKDICZPcH7KiyqNfjpSnZRMn0EgBqt9Wql+BAyTQPrCRsX6KKuSxmOfvzX77BDyDrexYKGaJRdz2hiIiISG+pWbzIKBJeHu52mwnHOASKAsQWaw7iQHBnhlrS2b9vcQ27Gcf1XIkbCBockgS864lVOooERUREpNc0HVRkFEkXkemOlE0Rr48TWRPp30EJkJ4ZatIzQ/ke32Qez/M5fkFmWige1jUeTKykTFNDRUREpF8oCBQZYdLVRWuX1OL1eLt1TsXaCubcOEfrBgdANOoWgDEGAqUp3jetig+zkjO4DwBDirmsI0ycCJHBHayIiIiMSAWDPQAR6R/+aX6akk3dPr5mWw0ly0pUUbSf+f3ZFUC9wFOwronxgY/yFlNYzQe5nGVuMRmlAkVERKQfKBMoMoIFfUEc4/4xNxi8Hi8Gw4xxMzo8J14fJ3hzUFnBAeSfW8gp8Ru4f9ynOIqnWMJS1sZ9ah4oIiIi/UJBoMgIFl0UJVAUwGAI+oJUXlJJqjzFlq9u6bS9RGOykbIVZQM4UiEYZNLr1fx1zMfwk6CMP/Fk5SS3qqiIiIhIH1IQKDKCpdcHpspTxBbHcqZ4RhdFKZ1e2uG58fo4xdcVKyM4kG64gRl7XuVhTqSILZxi/8Lz8QJ3AaGmhoqIiEgfUYsIESHRkCBwU4CmVO4awnQGUS0kBlA4DJWVrLezeD//ZI8znrUXryB486WDPTIREREZRtQiQkQ65Z/mp+rSKnzjfTnbLZZ4fZzpP5jOnBvn4FQ4hJeHlR3sT9EoBIPMZgN/d07CSTVxwvKPUTvnFK0PFBERkT6hIFBEADcQrLuyjpAvlCkmk9awp4GabTVYLFVbqrResB8l8BMmhmMslzo3cQ9nsYtxHFrzBx754NWDPTwREREZARQEikiO6KJop9VD0w3mtV6wf5SVQVUVWAsPN3+QxSznXs6kkGZO2ngrr5n9VDVURERE9oqCQBHJkZ0R7MyWnVuUEexjkQjE45BKpbcYYszjSq7n93wcg+UEVlFfuUVVQ0VERKTXFASKSF7RRVG8Hm+H+9PrBbc3bh/AUY1skYjbI97J/JfZApZnOYIIEX7Px0ng52T7V96M/1cVQ0VERKRXVB1URDqUaEhQtqKMyvpKgr4gjclGarfVYmn970bR+CLGeMbw2juv4Z/utqBY37CeuUVziS6K5rSlkK4lEm6Sr7ISSkrcbTU1UFq4gWOaHuUs7uET3MuhY6pZWx9m0qTBHa+IiIgMTZ1VB1UQKCLdlmhIELzZDQa7K+QLKRjcS5EI3HtnI1U1Ho7n71zEz1nECg48AF6sLGDixMEeoYiIiAw1ahEhIn3CP81PU7Kp6wOzxOvjlCwrUWuJXopEoKIC4jVeUnh4mJO4g/O5g/PYtBn2K27mpZcGe5QiIiIynCgIFJEeCfqCGEyPz1Mw2DuRiFspNHut4ErKuI8zuZ3z2bnLcNy7m9m5c1CHKSIiIsOIgkAR6ZHooihBX7DX51fWV6qqaC9EozAj07nD8gc+zv2cwa+5gB07DZMnNFFSuJE5BzfiOOoiISIiIh1TECgiPeKf5ie2OEbtklpCvhAe48E33tft89NVRdVnsGf8fqirg/JyCHlrcUhyH2cSpYxf8llSeNjQvD81Gwux1u01qC4SIiIiko+CQBHplXQw2Hx1c6avoGPc/6Q4xiHkCxFbHOswQNy6a6sygr0QIUK08SQCVOGhmd9zFn/hI/ycC7E4eEgCbq/BeNyqi4SIiIi0o+qgItIn2raTyK4I+sgrj3DCHSewJ7mn3XnTx05n+vjpJLYlKPAUkEwl1V6iK+Gwm+pLpQjzElUEOJPfcyIPcxG/wKEZiyHoXU9sT+lgj1ZEREQGgVpEiMiQcPBPDmbjWxu7PM4xDoGiALHFsQEY1TCUbiYYj5OY9QHKNi1nXaqUy1hGgHVczM+YxFsUU0eCUoIhQzTqTikVERGR0UEtIkRkSFj96dUUjSvq8riUTRGvjxNZE+n/QQ1Hfj/EYlBejn/DamLJAN9a8Ag/Dv2SL5ifcwsX8w5T2MDBWIzWB4qIiEgOZQJFZMCFl4eprK/Ekv+/P8oE9lIiQd2xZcyor+QWvshl3IxDMykKALeojNYIioiIjA7KBIrIkNJVm4lAUYDoougAjmiE8PsprovhKb+aTfsczU1cQooCPDQR8DcqABQRERFAQaCIDAL/ND9nhc7Ku2+MZwxnBs9kv4n7DfCoRpBIhC9MWsEmDuAWLiZJIcWvPceuXYM9MBERERkKFASKyKCILIy0aysxZ/ocPh76ONf88xrCy8OsfHnlII9yGIpEwBj8NQ/xA77JxfyMX/JZ/rXrKD4aqmbnzsEeoIiIiAw2BYEiMmiii6IEigIYDIGiAH8792/cecadrD5/NWMLxlK2oowT7jiBOTfOwalwCC8Pq8F8VyIRsBZCIXDc/8R/1tzGbVzA6g1+PnJyEzt2aG2giIjIaKbCMCIyJDUmG1n2xDK++vBX2xWQCflC6iPYlaw2EoRCcNll/G7xI5xnb2fCuBTv7CokFKK1dUQkoshQRERkBFFhGBEZdja/vZn/e/T/8lYQjdfHKVlWosxgZ7LaSBCNwo03co69k99wHjt3GQpoorISyk5upK44DBUVbhP6hD5PERGRkW5YBIHGmNONMT83xvzRGHPSYI9HRPpPoiFBeHmYkmUl1O+s7/TYyvpKylaoAV6nIhEoKyNVWQXAbzmX33AeAB7bRLymkJn1zxPmJRKVe9o1FFRyUEREZOTp9yDQGPMrY0ydMealNttPMcasM8bUGGO+3tk1rLUPWGsvAi4Azu7H4YrIICtbUUZlfWW3jrVY4vVxCq4p0JrBjkQiJOK7mGdfoIAm/s4H+RWf4U4+BUABTSQpoIoAZfaP7vTRSIREwk0MKkEoIiIy8vT7mkBjzPuB7cAd1tpDW7Z5gJeBE4HNwH+ARYAH+H6bS3zWWlvXct6PgDuttc90dk+tCRQZniJrIlSsrWi33WAomV6C1+MlXh/v8Hw1mc8vPKaGqsbZpPBgSFJIM8fyOJdyI5/idwA04cVDM6/6Dqe4LkY4DFVVkEq59WUCAXd2qYiIiAwPg7om0Fr7T2Bbm83HADXW2oS1thG4CzjNWvuitfbUNo864/oh8NeuAkARGb7ato1IC/qCPHjug8QWx6hdUotvvC/v+SmbIl4fz2QF12xYQ3h5mIJrCkZllrClWwTxxhJSeACweGjEi5k1iz8e/L/cxScxWArZQ0nhK6xYFHXPibsBILhfWxKEIiIiMgIM1prA/YFNWa83t2zryGXACcCZxpiL8x1gjPm8MeYpY8xT9fWdryMSkaEru21EyBeidkktscWxTCVQ/zQ/dVfW5Q0W09LTRD9w+weI18dJ2iRVW6pG3frB1m4RJt0tAsdxX6/ZOJtr/v5e7j34q9zF2XiwmH2KWfS3T2MxhLw1OI7NOkdBoIiIyEgxWEGgybOtw3mp1tpl1tojrbUXW2t/2sExt1prj7LWHuXz5c8SiMjQ55/mJ7Y4Rqo8lRP8tZUOFgG8Hm+X101nCYuvKx516wejUXc6J0DA30i08WQA/GVhVvy9mLrPX81fJn+STZsNC1/+Ga+yH9GmUwgU1LrnBNxrZCgaFBERGdYGKwjcDByY9foA4NVBGouIDEPpYNGWW/Z8aw+23BIsCnZ4vMHg9XjZumsrFktlfSXBm4OjYqpodreImPdd+BOr3B1VVVBWxhf2+xMfePuP/JUPsYkD+R/+hWObiTXOIYlDjDB+Eiy9XNViRERERoLBCgL/A8wxxsw2xniBTwJ/GqSxiMgIsfKclYR8obz7HOPQmGwkZd2FbhZLY7JxVE0VjRDJv9ivooLEkp/wvcJyfsO5bGM67+VRXmYODpZE5R7CwRRfXnYQ4fg9JJidCSBFRERk+BmIFhErgMeAucaYzcaYz1lrm4FLgQeBSuAea63qzonIXklnB2uX1GaCwUKnEICkTXZ4XnqqaGRNZCCGOXgiEXdxX+4CQbCWslWX83Dyg1zGTdzGp2mkkPfzT17kUMrsH6lqnN3aSoJoawBpjKaHioiIDDP93iJiMKhFhIikffPv3+T7j7TtPAMe48FiSdnU6GotkUi4Gbx43A0Ao1Eid/ipyOrMUcwb/JwL+SI/ZTdjeZMppCjI7PfQTKMZgxNU3wgREZGhalBbRIiIDKbvHf89Qr4QJqse1YTCCaw6f1WmsEyBU0BlfeWIXxsItFkgGAO/v12CcIuzDz+YfSv/8p7AFN6igCTGtFQKJclc1rkBYE61GBERERkuFASKyIgXXRQl6HOLxuw7cV8Mhg/f+WE+ffinCRYFaU41Z9pKlCwrGR3BYJspnDkVRAPw21X74a9/gn/yfmazngKa3X2+rUSXrHIDyDvuGOBBi4iISF/QdFARGXU2v72ZS/9yKX9c98e8+w2GoC84OqaHthGJtIkPr7qKNx6t4aS13+QlDuUBTqfMRt3jKozbiFBERESGnM6mgyoIFJFRyVrL/VX3c/bvz6Y51Zz3mPIF5UQWRgZ2YEPRnj1sO3kRp6z9Os/yLn67z1d41xt/4RBqMusK8efv5ygiIiKDQ2sCRUTaMMZwRvAM/nPRf5g2dlrOPsc4hHwhBYBpY8Ywva6KVZzAe3iMRW/8mACVhHmJmvgetYoQEREZZhQEisioNn/f+Wz72jbuPvNuvB5vZnu8Pt7p2sAR304iLRJx20BUVjKZd/gCP+VEHsZSQCUBTuOPahUhIiIyzCgIFBEBPhH+BG9//W3GF47PNJSvrK/k1N+dmnNcoiFBeHmYirUVo6eAjLUQCpEyDp/mdk7gYU7jASweKglQSjWXL7HtFxQqKBQRERmSFASKiLT4/iPfZ2fTzsxri6VySyWX/eWyzLayFWVUbakCoGpLFWUrRslUyGgUJxhgLi/zdX7A4TzHIn6HxcNGDuLhW2rcHoTZDQezn4uIiMiQocIwIiJZwsvDVG2pymQD08YXjs8JELONlgIykUhuXHcZS9nJBH7JhXjZzZtMoZAk23xzKX48CiUlqh4qIiIySFQYRkSkm6KLopkm8l6PN9NkfmfTTjzGw+yps3GM+5/O0VZAJj0z1JZHALiZS9nNGC5lGY2M5VKWY7AU18fhkEPck8JhN0MoIiIiQ4aCQBGRLP5pfmKLYyyYtYDGZCOW1kxW0iZZ/+Z6Jo+ZDECBU0BlfeXoWBuYLRJhyRIIeKpZwSI2sz9f53v8is/xUf7EHryQTLrHVlaqeqiIiMgQoyBQRCSPNResIeQLZbJ+2d7c/SZej5emZBMWO7rWBrZYuhRiL3uZzQb+xOk8yvu4lm/wFz7CR/kTOxjvHmitqoeKiIgMMQoCRUQ6EF0UZca4GXn3ZWcJUzZFvD5O8XXFoyYjGImAKfFTSykpPPyL93MPZ3MDX2IVJ3ASD/EmU9yDQyE3GFQQKCIiMiQoCBQR6YB/mp+6K+tyMoLpdYDpdYPZtu7aOmoyglmdI3Bwp36+yDwe4mTu5mz+w9F8gH/wBsWtmcCFCxUHioiIDAEKAkVEupAuFmMwBIoCRBdF+fM5f8Y33pdzXDojaCrMqGkmH41CIOQBwF/qoaYwyP9yA7dwMeuYS4B1OKSYU7ieOWtv5bsVzYTnNKpWjIiIyCBSECgi0oV0sZhUeYrY4hj+af6cLGG6gmi2e2L3kGhIZILBkRoU+v0Qi7nPvV5IJGexmYP4Gj/kx/wvFhjPDmqaZlHDHJIUUFXjUa0YERGRQaQ+gSIieyHRkKBsRRnx+jgFTgHNqebMvkKnkKZUE16Pl8ZkIyFfiOiiKP5p/kEccd9q2zsw23h28EO+yne4mu1MZCcT2h1TvmANkTUL+3WMIiIio1FnfQIVBIqI9IGFty1k7ca1nR7jGIdAUYDY4tgAjWpghcNQVQWpFIAFDIXs4Qd8nZ/wv9ThYw9jcUgRCHkyGUQRERHpe2oWLyLSz9ItJfJNDU1LrxnsamrocJw6Gom49V/cABBo+RyaGMOVXMcVXM8sXsFLIwdN3EY0OlgjFREREQWBIiJ9JLooStAXxGDwerztAkKDIVgUJLIwkvf8REOC8PIwFWsrhl0D+nS1UFseAVorhjokCbCOL3Ej/+J/CBHnv9un8PT/fIl0dRhVDBURERlYCgJFRPpIdgGZyksqCfqCgLs2EMBiGVswlodqHyK8PIxT4eQEe2UryqjaUgUwfBvQRyIsWQIB73o8NBOgiqg5DUIhikM+/sEHeTdP8MlXf8QP3n0/4eI6KircqaSqGCoiIjIwtCZQRGQAWGtZ8dIKlvx1Cdt2bXO3YTEYCj2FNCYb8563YNYC1lywZgBH2kcSCeqOLaO4Pt5u107G8XHu4298iIm8zXYm45gUgaCjdYIiIiJ9RGsCRUQGmTGGc+adw2fmfwbb8j9wA8F8AWB6SunajWuH3dRQAPx+iutiUF7eMk/U7SyfMg7j2cUDnMYZ3Md2JjOZt0hZJ9NTXtNDRURE+peCQBGRAXTdSdd1WUAGoNBTmGk3MWynhkJuRHfCCTjWrRzjpYmx7OJ8budtpjCNrQSJYTFEiOS9lIiIiPQNBYEiIgPshNknZNYLZgeDBoNvvA+AxmQjqZaAqbtVRYe8pUtJ1FrCvjcopIknOYatTOOLLKeBGQRYR1PwMDj/fHeRoONosaCIiEg/UBAoIjJA0tU/lz25DIAlxyyhZkkNB0w+IHPMpw//NE3fbiLkC+EY9z/RjnEI+UKZqqLDORgsK4OqLUUkKSBBCesp4WYu4Vq+yf2cwYcqrydYspuC+POE7YskKvfAsccO9rBFRERGFAWBIiIDpG31z1XrV+Gf5mfT/27iK+/5CmeFz+L6x67nPb98D9efdD2BogAAgaIA0UXRYdtCIhJx1/oZ09JL0Lp/9aTwsI65pHD4Jt9nORezmuOpp5gkBVQRoMz+EerrWy9w+eWD+2ZERERGAFUHFRHpZ5E1ESrWVuTdV76gPKdv4O/jv2fxnxfz1p63KF9Qzu7m3VzzgWsAKL6umK27tpKyKRzjECgKEFs8vMpphsNQVeU2lXccCPgb+cdb72J6/TqacbieK7mGcibxDtuYgYdmmijEACnj4AQDqISoiIhI11QdVERkEEUWRrDlttMpnmlnhs4ktjjGxwIf46rVV/GX6r/w1+q/UnxdMfU769utEzQVZlhND41GIeAmOAkEIPrhW9hev4vDeZ6J7ORHXEE5EfYwhqk0MJtEZtWkY1NkSoiqjKiIiEivKRMoIjJAEg0JylaUEa+PE/KFiC6K4p/m7/D4++L3sfgvi6nbUYfBZNpKAJlM4Fmhs9oFksNBJNIaw2VnBw1JxrGTr/NDfsL/0kghxdTxMCfhdza6kaMygSIiIl1SJlBEZAjwT/MTWxyjfEE5scWxTgNAgBfrXqRuRx1ATgAIMHXMVBqTjcNufWC2nHWCboITi4c9jOMOzudiljOVt3id/TiJB9ky1Q+NjaoaKiIispeUCRQRGeLCy8NU1ldmAsHxheOZNWUW67auG9brA9OyM4G5LJexjIc4mfXMxkMTB7ORlZQpKygiItIFZQJFRIax6KJopq/g5DGT2dm0k8otlcN+fWBaep2gMRAKQW2t+9UhxY1czpH8h8N4gUbGsokDKaGWcOoFEvFdrelEY5QdFBER6SZlAkVEhglTYbDllgeqHuDMe84kaZPudgxBX3DYZgLbiUS4vGIaD3MCLzOXJB4+xF/YxXjW8AGmsZW3mEqAKl7i0EzhGLfcqLKDIiIioEygiMiIUL6gHIDTA6fzxIVPMHnMZADGFIzhRyf9aDCH1rciEZbay4mHPsEexhAizoOcwm7GcCpRGphBEfVUMbc1AAR3Pml29VD1FBQREclLQaCIyDCRXQX0yJlH8tbX3+IT4U8wyTuJ0+46jRN/cyLJVHLwBtjXolE8oQBRygh41/M476GGEhZxJ3XsyxxqSHb015jjwKpVAzteERGRYUJBoIjIMHb3mXfz53P+zBjPGFYlVjHp+5P4+/q/D/aw+obfD7EYftYT21NKba3BCYVYy/v5Aj9lHQE+yV3swdv+3LZZQfUUFBERyVAQKCIyzF3wxwvY3rgdgF3NuzjxjhO56cmbMoVjRoqWmJCLyg/k68+czZecpfyeswh7q3nxBdtSTabNX2vpSjMKAkVERDIUBIqIDGORNRHi9fGcPoIWy2V/vYzSZaW88tYrgzi6/hGJwEfOncat9vMs4SdsaJzJR9+1iS23rXQLw2SrqoKyskEZp4iIyFClIFBEZBiLLIwQ8oVwjPufc8c4BIuClB1SRv3Oeg5dfii/evZXZFeC7k4biaHYaiISyW0wv9OO42Yu5Qv8jNeTRRx2zBgK4s+Ryi4Xo2mhIiIi7SgIFBEZ5qKLogSK3AyYf5qfplQT0Zej7DtxX0K+EJ/70+coW1HG45sfJ7w8TMXaCsLLwyQa2vfUSzQkujxmsEQiYK37SM/8TFLALSzmk87d7GACM0wDqw84v3VaqOO4B6dPVBAoIiKiPoEiIiNFZE2Ee+P3UrWlipRN4RiHuTPmMmvKLNZsXENjshFrLRaLYxwCRYF2vQXDy8M55+c7ZsAZ4wZwWRIJd5ZnPA5eLzQ3w0dTf+DfvJdGxvDQQRdy9KY/QDDodqP3+wdp8CIiIoNDfQJFREaJeH08UxAmZVNUbqnkb7V/Y3fzblI2lVk7mLIp4vXxnGmf6fWF2efH6+OYCjPkpofecYcbAAI0NrqzPh/gDELEmEoDH3jlNmbaTUTOiikAFBERaUNBoIjICNF2fWA2xzgUOoWYrPVyB04+MKf3YL71hSFfCFs+yDNGysvbbco3NRRgDR9kGtuYxUa2UsS88TXujkQCwmH3wHDYfS0iIjJKKQgUERlBstcHZkvZFE2ppkwmcGzBWDa9vYn5P53Pm7vfzHt+oChAdFEUgIq1FTnXG9DMYBfr+KLR3KKgz3IkTXg5zLzIJ742m/d7HqGg5CDC8XtI2INVMVREREY9BYEiIiOIf5o/s4avs6zeW19/i2+//9s8/8bzzLtlHg/VPtTu/PTX8PJw5mu6UEzboHCwRCJQUtI6NTStmjlssvvzXh7lX6n3cTjPEydAkDgFqT1uQGj8qhgqIiKjkoJAEZERqHxBeYdZvfIF5Wx+ezP3Vd4HQN2OOk7+7cl8ceUXM03nyxe4UzDLVpRRtaUKgKotVZStGFoZtHbTQk0qs+8N9uMFDuWDrOIZjuRonqIRL0kKqCJACQki5aoYKiIio4+qg4qIjHCmwrRb15ddBdRgmDZuGtt2bcM/zc+vT/s1q9ev7jDb5xvvo35nPSFfiBNmn8DSDy0F3Cmi2WsMB1p2xdDSUndbTY2lgEaOZzUP8iGO5kme4kgsnpxzy8sVC4qIyMii6qAiIqNYOqsHbqBmKkxOFVCLZduubQBsemsTC29byDt73mHXVbvaTSn1erxs2bkFcDODy55clrn2YE8R9fshFnMDOq83XfvF0IyXBzmJD/Nn/sMxHM7zjGE3Xq87GzQUgvPPH9Shi4iIDChlAkVERqHsTGA2g2Hq2Kk07G4g5Avxg+N/wNf//nXi9fEOrpQr5AsRXRTFP23g2zJEIlDRRRx6PA/zDz5IgCpeZX/eZCqO4xaWiQ1yO0QREZG+pEygiIjkiC6KMmPcjHbbLZaG3Q18at6naNjVwMfv+TifmvcpAGqX1OL1eDu97mCuG+yobYTjgM/nPv87J3IKf6WWEqazlX15lVQK4nFLzIRpNgXUFauFhIiIjGwKAkVERiH/ND91V9YB+auI/vaM3/LiF1/kY8GPcdXqqzhg8gGc/NuTaUo2dXrdodJgPrttRCAAjz/eGiD+OfQ1PsBqtlCEAUp5GYDTuZ8gcWbWP084mFIcKCIiI5aCQBGRUayzKqIzxs/gro/fxRmBM9j89mZqttVk+gx2ZKg0mM9eHxg7K4K/xLgLAI2BeJybuZQjeYpmCthKEfN5lhrmUMMct3po42zKSmJqISEiIiOS1gSKiAiQv4po2ua3NxO4KcCOph2ZbV6Pl8ZkI6XT3VKcNdtqctYEZl9vsCuHZkQi1FUsZwZb8ZDiQDZggK0U8S6e4VH+J+9pqh4qIiLDjdYEiojIXjlg8gE8d/Fz7Dtx38y2H5/8YwCqL6um+rJqoOMG84NdOTQjEqG49nE8ITfzOZnt7GA8B7CZJziWE3gQAIckIWLYUBhbm1AAKCIiI4qCQBERAXJbSeRTOr2U1778GgD7TdyPS/5yCQBv7X4r57hjf3FspsF8ZX0lwZuDQGtAOOjSc0WXLCFaegX7UsdmZhJ2KlnFyZxKlACVRCmDqiq3+aCIiMgIoiBQREQAuj1ds3xBORu/tJGr3381BsO8W+axev1qFsxagKkw1O+sz+lB2JhsBCBeH+fYXxzbX8PvuaVL8Vc/SKz8XnYwmcdTx3AG97GSMj7E3zidP1CQ2kM4fg8J49f6QBERGTG0JlBERHrtyf8+yXn3n8fLW1/mM/M/w2ObH8tkAbuyYNYC1lywpn8H2EPJZsslk+7gZ7s/zSn8lX/xPnYxnkDIoz6CIiIyrGhNoIiI9Itj9j+GZ7/wLJccfQm/fu7X7QLA7L6CjnHwjfdRu6SWkC/E2o1ruz1FtL/bTUQibqKvoNDw4O7/YT7P8Dc+xNE8hY864vHW4qJKBoqIyHCnTKCIiPRaZE2k20Vf0tVEvR4vzalmUjaFYxwCRYFMQZmOdFa5tD+Ew0A8TiUBjuY/NDCVwtBcolF3SaGIiMhQp0ygiIj0i8jCCLbcYsstIV8Ig8nsK3AKMlVDgczawMZkY2bN4FBpLt9WNAqEQsxiI8/yLrw0MSn+b4IluykwzYTnNKqZvIiIDFsKAkVEpE9EF0UJ+txKoPtP2h+P8TDnxjldnucb78OWWyILI0MiEIxEoKQEdsUTbOJAprGN9czmDfbjA6x2m8nXeCgp0dRQEREZnhQEiohIn/BP82emdW6+YjN/P//vTCicAIBp+V/6eXqtYMgX4vELHyfRkCC8PEzF2opBbyURiYC1kPAdy1zWsQUfBTTRwFSe5UgWcScpPHho5qyKMH6T0FpBEREZVhQEiohIv/j8ys+zs2kn4LaKsLhr+oK+IJWXVAJuc3n/ND9lK8oyRWWqtlRRtsLtzZcODmGA+gymK8QYA/X1RCkjQBU7mEgzDhbDSsr4IjcToJKwU0UiVIa1CgJFRGT4UBAoIiJ9Kt0vMF4fzwR+2Y7Y7wimj5ueOS59bL51gtmN57ODw36TTgNaC6EQfmcjMQ6lmUJe4F1MZwsT2M4vuZCv80NIpcgpHaqUoIiIDAMKAkVEpE+tuWBNplCMY9y/ZhzjECwKsmDWAla8uILDbjmM8gXl2HJL+YLydsemZTeeH/AiMtEoBAJuYBcK4a9dRZUNETvk4xzFU5zPHfyUL7Qe7zgQCikIFBGRIU9BoIiI9IvooiiBogAAgaIAK89ZyZoL1vDvz/2bcYXjOP6O47nyoSupWFuRc2yBUwC46wVLp5fmBIchXyhTRKbf+f0Qi7nZvljMfR2JMP3lx7mVi5hNgi/yU77EDWzkgPZZQQWDIiIyRKlPoIiI9Kt8Pf52NO7g8ys/z+9e/B3gFpV5+LyHOfYXx7J119ZMD0H/ND9ej5d4fZyQL0R0URT/tMFv1BcOQ2UlFNvXeYN9uYBf8TyHsYfxRCnDHxqLmgqKiMhgUp9AEREZUiZ4J/Dc689lKoYmGhKULCtpN/2zZlsN8fo4AGeFzhq0ADC7XowxbsLPWniDfZnBFm7js8xmI/uziTL+RCpeCcGgO0U0HEZNBUVEZChRECgiIgMmsiaSUwwmX+GYtOzpn0DOFNCB7ieYXS+mvDx331aKmMxb/IGPk8LDu3gaB4ttbHRPiMfh2GMHdLwiIiKdURAoIiIDJrIwgi237QrHGAz7TtyXsQVj8RgP4K4jjC6K5gR8nfUTHKjAMBKB2lq3Bowx4PXCdjOJseziH3yAdQSpp6glx9mivl5rBUVEZMhQECgiIoMiuxhM0Bfk0c8+yrNfeJbD9z0cgEN9h/KR332EirUVAFz+18vz9hMcjEbz2TVjKishEHTYzVhmef5LjDDv4xE2MMs9OF01NJ1KVBAoIiKDTIVhRESkX0XWRDqt5tm2cExjspHImgjff+T7Pb6XYxwCRQFii2O9GepeiUTcx78ftXzkg7sY39jAd/gWX+InHFg6juiDXtWJERGRAdNZYRgFgSIiMqiyg8TImkgm89cXyheUD0w7iRaJBJSVucsAJ7CDQhr5IV/lanMtM4LFxAY+NhURkVFK1UFFRGTIyin40rJmsHxBeccnZFlyzJJ2jebTxWQGrJ9gi0gESkrcABBgBxNoopDLWcY19iqS8UpiJozfJDQjVEREBpWCQBERGXIiCyPULqkl5Au127ffxP24+v1XY8stSz+0lBNmn5DTlD66KDpw44y01nupyJPA3MFEPCRZzC1cxo2EiZPwBolco9YRIiIyeBQEiojIkOSf5s+s7UsHg+MLx/Pa9td4qf4lnn71acLLwyx7chngZgVji2MD2kswu3WEtW79F8e4fQ4NSbzsZhdjmUoDl7Kc6/gKZLeOKClRoRgRERlwCgJFRGRIK19QngkG3/nGO/zfCf9HdF2UY395LJX1lYBbKXTV+lWDOUwAolG3UqgxEAx5qFzyM5IU8gjvo4g6vsp1XMH1XE2EApoI8xKJijvUPkJERAaUCsOIiMiwsPC2hazduLZbxw50QZhOJRKEgykqG2fjIUUzhZzP7czhZSJEmMvLxEKfcCNIlQ8VEZE+osIwIiIy7K25YE2m4Istt5l1gGn+af7MvsGWvVbQlPiJN5Zi8dBMIYYkd/BpHuc9/IwvkGA2yXglu0uCNJsC6oq1VlBERPqXgkARERmW/nzOnzNrBQucAja/tZlv/v2bhG4OZRrHX/7XyzPHR9ZEBmxsna0VtDhAir/wYX7JhTzA6XiwjKWRApIU12utoIiI9C9NBxURkWEtsibCJUdfwudXfp4Hqh5ot790eikANdtqCPlCRBdFB7R4DLT2D6ysdOM7tm0lsW0qBsshrGMZS/gB3+BWPo+f9e0vUF6uoFBERHqks+mgBQM9GBERkb5WfH1xh/tqttVknldtqaJsRVmm0MxA8ftp0yh+BiQSHDx3DC83z+Uz/Jo7OZcr+T/uc86GQKDtCSIiIn1G00FFRGRY60mD+ZRNEa+PYyoMpsIM6BRRaL9WcGPz/iQp4DVmcjoPcCk3cWrqj3jjz6pYqIiI9BtNBxURkREj0ZCgbEUZ8fp4p8cN1rTQtsJhqKqCVAo8NOOlkfv4OA+YM3gkeJGSgSIi0muqDioiIqNCvgbzBU77lQ/paaGDKRJx+8Wn3HoxJClgD2P4KH9ioV3NufFvEDMh/CahrKCIiPQpBYEiIjLiZDeYb041t9vfdlroYE0NbVtB1BoPzXg4hxVMYAdhKkkEPoK1CgJFRKTvKAgUEZERJ90ovnxBeaZ34MxJM3OO8Xq8gJsxrF1SO3jN5VsWCkbjfoI2hkOKCbzD5Szj21zDE1WT2GBmtS4mVDQoIiJ7SWsCRURkVEg0JDj5tyfnVAsFcIxDoCgw4BVDOxQOE4rfQyVBwOFz/JyvFi7lkOo/w6xZgz06EREZJrQmUERERr07nr+jXQAIXVcMTTQkCC8P41Q4hJeHSTQk+nxs2VVD/fEoLzOX9F/Rv+Qivtz0fXYdPJc9ppBmU0BdcdhtPigiItILCgJFRGRUSLeSsOWWkC+EwWT2eYwHcKeGNuxqyDmvbEUZVVuqsNh+KyiTvT4wYf3MDRXgtPwNbUixkjI+yGq24OMylrF9yx63+7yIiEgvKAgUEZFRJ7ooStAXBNzpoEmbBCBeH2fZk8tyCsbE6+OkrFvCs9/7DGatDwykYhhSBKlkH17jcd7DB1nNt/guv7HnuKVFM00HtVZQRES6T2sCRURkVImsiVCxtqJbx5YvKOfe+L1UbakiZVMDsn4wEoGKrOG9RJjDeZ4kBYBlX17nEd7HLDYAsI4Aq5ZEuXzp4PY8FBGRoaWzNYEKAkVEZNTqSUAIg9RkPpEgHExR1TibFB7AMpm3WcsCSqlhotkJwSDqLC8iItlUGEZERCSPyMIItUtqM43l80m3kEj3HhzQABDA7yd68Z8JUAW4/3D7NpN5N4/zFEe5CwnbTg3V9FAREemEMoEiIiLA5X+9nFXrVxGvj1PgFGSazBsMgaIA8UviAz6mtlND23JIcg9ncQb3Y4AkDlt9AYrrlBUUERntlAkUERHpwrRx04jXu4FeOgAEsFgqt1RmisEUX1eMqTCM+e6Yfm0bAblVQ62FUAgc4xapMSQZwx4+wb38ggsB8JCiuF5ZQRER6ZwygSIiIm2El4czxWAMBmMMk8dMZqJ3Iq++82qmWigMbLP5RMLtDFFZCSUlkEpBIpECHCr4Nt/mu5gDD4RXXun3sYiIyNCmTKCIiEgPRBdFCRQFADcTmLIp3tz9Jpvf3pwTAEL7thH90jqihd/v1n9JpcDrhQ0bwP2rPEU53+Fyz82kNr8KN9/cL/cXEZGRQZlAERGRDkTWRIgsjADQlGxi5g0z2bJzS84xA5UJ7Gp9YNr/8E8e5kR2MJ7mov0ofmKlGz2KiMiookygiIhIL6QDQIBCTyFPXPgEs6bMyjlm7oy5RBdF+38skTzrA1v+FnccKCx0l//9i/czn+copJm1W8KES3Yyx1RTYJoJz2kk0T/LF0VEZBhRECgiItJN/ml+NnxpA19779c4O3w2AJPHTOaU356CU+Ew58Y5zLlxzoAUjInH3Wmh4H5tanKDQ4AqApRSzfv5Jz/jYrYxjSQFVNV4KCuJQTiMokERkdFL00FFRER6wVrLHc/fwWf/9Nl26wRhsKeJWsAwgXd4kndTQDMf4c/UMAcPzexhDFUEOJQY5eUqHioiMhJpOqiIiEgfiqyJ4FzjcMEfL8gbAMLAFYxJTxOtrU23kLB4acSQZAeTmM9zbOYAHudY3s8a5rIODynCxLEYIhVqIyEiMtooEygiIrIXwsvDVNZXYsn9+9RgKPQU0phsJOQLEV0UxT9tYAq0ZFpJxC0ekljgF3yOc1jBm0yl2NkKgYBbalREREYkZQJFRET6SXRRlKAvCMDMSTPxGA/gTgdtSjYBULWlirIVZQM2pkwrCWvYsLmA4KEFfL7g19xTeC7F1BNJfZvwy/eTMCVaHygiMgopEygiItIHFt62kLUb13b7+PIF5TnVR/tTQwOcdhr861+Wc/ktv+F8fsO5/IgreM45SllBEZERqLNMYMFAD0ZERGQkWnPBmsxzay0zb5jJ69tfz2zzerw0JZsI+oIDMjU0f8EYw285jwSzeYT/YRYbIZWiOb6OQgPlRIjQcpIqxoiIjFjKBIqIiPSDREOCE39zYrs2EQNVNTSfUAgqK93nJdTwAvMYz24oLYXq6gEfj4iI9B+tCRQRERlAkTURSpaV5O0TOFBVQ/NZuRKCQQBLLaWczh/Z6UyErVvh0UcHZAwiIjL4lAkUERHpZx1VEB3MrOB3vgNXXw1j2cVjvIfDeR6z//7wz3+6lWVERGRYUyZQRERkkETWRIjXx9sFgDC4WcG77nLbA+5mHEfxFH/mw/Df/8Kxx7qNB0VEZMRSJlBERGQAhW4OUbWlKhMUzp46m8Tl/duiIX+RmFyGFDHCBKmigSlM4h22MYNi6lUkRkRkGNrrTKAx5lJjzLS+HZaIiMjos/KclZm+ggVOAZve3sR1j15Hyqb67Z6RiJvcy36EQuC0/BZgSOIhyft4hMec9zKNtyggRbHZ4h6oAFBEZETp7nTQfYH/GGPuMcacYowx/TkoERGRkco/zU9scQxbbqn7Sh2nzT2Nr676Kqf89hSe2PwE4eVhnAqH8PJw3sIyfSUahcCMOhySBKliFccznW0cn3qIlXzEPchaiMfdeaPphwJCEZFhr9vTQVsCv5OAzwBHAfcAv7TW1vbf8HpH00FFRGQ4SDQkOPV3p1K5pRKDwTEOKZvCYjEYCj2FA9pbsK4OPjLrJZ7dHeDnXMRnuA08HrdYzHHH9eu9RUSkb/VJYRjrRouvtzyagWnA740x/9cnoxQRERlljv3FsVRucRv3WSxJm8ysFbRYGpONWCzx+jgly0r6rWhMIgHhMOy7L7wz8xDeM+55Psuv+cGU72ONA+99LxxwAKxZ4x5YUOB+TfTvWkYREekf3coEGmOWAJ8GtgC/AB6w1jYZYxyg2lpb0r/D7BllAkVEZCiKrIlQsbaLCi3dVL6gnMjCSJ9cKxyGqipIpdx1grNnQ309vP027Mur1FDKBHa5O9OLCh0HAgGIDXx7CxER6VpnmcDuBoHX4E793JhnX9BaW7n3w+w7CgJFRGQ4CC8PU7WlipRN4RgH/zQ/Xo+XeH085ziDIegL9lk/we5UC802lQZeZ1/G0Nj5gaoiKiIyZPTFdNDZbQNAY8xvAIZaACgiIjJcRBdFCRQFAAgUBXjw3AczRWMe/eyjjC8cD8CkMZO484w7AXcd4d4Wj8muFlpe3vXxbzKNseziGr6V2+3Q63WLxYRCUFurAFBEZJjobibwGWvtEVmvPcCL1tpQfw6ut5QJFBGR4SSyJpJ3amcyleT/Hv0/vv2Pb+MYh6ZUE16Pl+ZUcyZ7GCgK7HWGMJGAsjKorIRgEBob3W2pFJAJ+wxg+bTnt9yWPL9lk9HUUBGRIarXmUBjzDeMMe8Ahxlj3m55vAPUAX/sh7GKiIiMOh2t7fM4HvYk95C0SZpSTQA0JhszPQVTNkW8Po6pMJlHb4rH+P1u/Hb11W5HiJqadAAIbvBnMs9vT57HIVS5oWH6H5JTqfatJNROQkRkyOpuJvD71tpvDMB4+oQygSIiMpz1tIBMXxaJaSu7aIwxbtJvxgz469iPccQrD7gHGeNmAuPxTq8lIiIDZ28ygYGWp/caY45o++jzkYqIiAiRhRFsuc08Qr4Qjmn/V/aBkw+kdkltvwWA0NJUPuDGecEg/OpX0NAAR77yB0o863mHSW5GcOZM9wDHUfsIEZEhrtNMoDHmVmvt540x/8iz21prP9h/Q+s9ZQJFRGQkSTQkKFtRRmV9JUFfkGWnLOMrD3+F515/junjprNt1zZCvtCANJQPh921g+lfH2buZ/nvpd+Hq65qPUhrBEVEBt1et4gYbhQEiojISLeneQ8zb5jJtl3bgL5vIwE9ayVxKT9hKVfg0MnvFWohISIyYPqiT+DzwF3APdba2j4eX59TECgiIiNRb5rN93WGMHuNYK4UH+MP/IGz3Jf77QfTpsG6dTB3rjuv1N+/WUoREWnVF0HgLODslkcKuBs3IHylLwfaVxQEiojIaJDdbL4jfdVGontZwRQNTGUMTYxlNxbjZgY1PVREZMDtdbN4a+1Ga+3/WWuPBM4BDgPW9+EYRUREpAciayLE6+OdBoDQd20kshvMpx+hEDjGvb8hBRj24w1msIW/cTIOlj9SxsTUm4Tj95AwfrWQEBEZArq9JtAYczDwCdxsYBK421r7o/4bWu8pEygiIqNNR1nB/lgrmJZuMp+e8dnQAK+9lr5vM+VcQznf4T8cxcf5PZNCs5QMFBEZIHudCTTGPAH8AfAAZ1lrjxmqAaCIiMhoFF0UJVDkdnaaPXU24wvHAzDRO5Hffuy3fX6/SARKStzWgMmk+zUdAAJYCohQwWk8QIAqnuTdTI//k5gJ4zcJJQNFRAZRd9cEBqy1VQMwnj6hTKCIiIxWkTURIgsjpGyKnzz+E76+6uvsN2k/Vnx8BccdeFy/3jtf0RhjIGhj/JHTmMVGCk3SjR69XhWNERHpR70uDGOMOdda+1tjzBX59ltrb+ijMfYpBYEiIiKu//z3P3zyvk+yvmE9137wWr72vq/lbTzfFxIJKDu2jsr6IubwMgUkiRPCkOJYnmAVJzCeXR1fQC0kRET6zN5MB53Q8nVSnsfEPhuhiIiI9IsZ42dQ6BRisXxz9Td5/6/fz+vbX++Xe/n9EKsrJmUd1tkA/9kR5rS567B4+B/+xavsy6+4AIAnOYoNHJR7gYoKFY4RERkA3Z0O+l5r7aNdbRsqlAkUERFxZReMMRgAfBN83HH6HZxcenK/3z+ZhMs+/Ra33DmFybzJO0xiEXfxCy6kgWnMNK9DMKj2ESIifWyvC8MAN3Zzm4iIiAyiyJpITjuI7DYStuV/dTvqOOXOUzAVhm+v/navWkZ0OoZIazKvoABuuXMKAG8zFYuH3/Ep3scjJPGQsnBZ/AtK/omIDKCu1gS+BzgO+BLw46xdk4GPWWsP79fR9ZIygSIiIq6uGsqPKxjHruZdhHwhooui+Kf1X4GWAw6A//4363XhG9zV9HHey6M0TJvNtLc2uk3lVShGRGSv7U0m0Iu79q+A3PWAbwNn9uUgRUREpG91p6H8rma3UEu8Pk7JspI+zwpm++c/Ydas1tf/bd6HD7CaW7mIaQ3r3bKilZVu80EREek33V0TOMtau3EAxtMnlAkUERHJFVkToWJtRY/OKV9QTmRhZO/vHXFrvnTmi9zMci4FoBnDOoLMZR3bmE4x9S0DUvVQEZHu6nWLiKwLPIzbJP7NltfTgLustf2/orwXFASKiIi0l2hIULaijMr6SoK+II3JRhINiXaZQoMh6AsSW9x/xVrmzoWXX07fL0WASuIcmjUIA9aC47hTRFU4RkSkR/qiMExROgAEsNY2AMV9MDYREREZIP5pfmKLY1y94Gri9XFqttXknSpqscTr4zkFZvp6muhf/wqBorqWuznM4hVCxJjJZp7kaDcABHeKaDye2zpCVWRERPZKdzOBT+MWgnml5fUs4H5r7RH9PL5eUSZQRESkZ/IVkPnkoZ/kZ6f+jMljJvfbfa2Fa67JjecK2cNPuZjPcpu7Yc6c1rShiIh0S19kAq8CHjHG/MYY8xvgn8A3+mqAIiIiMriii6IEigLuVNCiIF9+z5e5N3YvR956JM++9mzmuERDgvDyME6FQ3h5mERDosf3ym4h4TjtE3pNjOHz/JwkDkkctlRvY6FZo+SfiEgf6VYmEMAYUwQcCxjgMWvtlv4c2N5QJlBERGTvPfLKI3zy95+kfkc908dP5/Xtr+P1eGlONZOyKRzjECgK9MnawXDYLQya/rXk4INh/XrcjR/7GNTUwA9/CFdc4UaPIiLSqb3OBBpjDHAKcIS1NgqMN8Yc04djFBERkSFmVWIV/33nvzSmGnl9++sANCYbM1NGUzbVbu1gb9cPRqMQDLrPvV549VVYsQJ34913w/jx8JWvwJQp8OyzkEi4kWNBgfs10fOMpIjIaFXQzeOWAyngg8A1wDvAfcDR/TQuERERGQSD1UrC728tALp1K5xxBpxzjpsA/NaKczDbt7s733kH3v1ut/P8xo1u4ZiqKre3oCqIioh0S3fXBL7bWnsJsBsy1UG9/TYqERERGRSRhRFql9QS8oUwGEK+EKXTS3FM+18ZgkVBapfU9kkvQffmETCGGUWGn/4zwBTe5OqrYWrlo8Tt3NbjmprcuaKpliI2HVUQ1QJCEZG8uhsENhljPIAFMMb4cDODA8IYM8EY87Qx5tSBuqeIiMhoVbaijKotVVgsVVuqADJFY+bOmMvxs4/PbJs2dlrf3TgScRcFWsuZoSreNlMBeJspHMkz/JeZbiWZuXNh4sTW84yB0lIIhcDjcb/W1ioIFBHpQHeDwGXA/UCxMeZa4BHge12dZIz5lTGmzhjzUpvtpxhj1hljaowxX+/G/b8G3NPNsYqIiEg3RdZE2q3pi9fHc9b91WyrIV4fx2JZt3Udf1//d04qOYnoy1GOuPUInvzvk3s3hkj7JF483lokBgy7Gcu7eIYDUhsw6yrxbt9CFYe4/zpdUAC7d7vTQpPJ1umhIiKSV0+qgwaA43Grg/7dWlvZjXPeD2wH7rDWHtqyzQO8DJwIbAb+AywCPMD321zis8BhQBEwFthirV3Z1X1VHVRERKT3snsGdlYB9InNT3D278/mv+/8lx+e8EP+99j/xfRR5c5w2I3lUik3+XfwwbBpEzQ1WcDgkCRAFTEO7dmFy8uVIRSRUaHX1UGNMdPTD6AOWAH8DnijZVunrLX/BLa12XwMUGOtTVhrG4G7gNOstS9aa09t86gDPoDbmuIc4CJj8ixKEBERkT6T3TMwUBQguiia97h3H/Bunv3Cs5x6yKl8+aEvc9pdp7FtV9u/9nsuEnEzgdlL/hIJdymg+2/RkMJDnDAGywy28A8WtF7AGLfEqNPyK4PjuFNErVUAKCJCF5lAY8x63HWA6X/WSx9sAGut9Xd5A2MOBlZmZQLPBE6x1l7Y8vo83MIzl3ZxnQvoJBNojPk88HmAgw466MiNGzd2NTQRERHZS4mGBKf+7lQqt7gThPabuB/3feI+3nPge/r8Xm17CRoD2BQWBw9N/ICv8xVu6NlFlRkUkRFqb/oEntcS6AWttbOttf6Wx+zuBIAdjSfPti7npFprb+tsKqi19lZr7VHW2qN8Pl8vhyYiIiI9UbaijHVb12Vev7b9NY771XHsc/0+1Gyr6ZN7pNcM5q4TbKkh0/KrTJJCruRHHMV/eJk5WFor2KVMViaw7UMBoIiMQl0FgUtbvv67D++5GTgw6/UBwKt9eH0RERHpB10VkclWt6OOOTfO6XXz+Jz7RvLHb7NmtR5jjBvnPWWP4pB3nsGcdVbmlxynxO92oxcREaDrILDJGPNr4ABjzLK2j17e8z/AHGPMbGOMF/gk8KdeXktEREQGSGRhBFtuseWW8gXl3T6vYm1FJmjc24Aw2+rVbmcIcIPCk05qWUc4cSLcc4/7mDoVXn8d/vWv3DSiiMgo1lUQeCrwILALeDrPo1PGmBXAY8BcY8xmY8znrLXNwKUt160E7rHWti85JiIiIkNWOiBMN5YH8Hq8mabyjnEocAoyxxsMJdNKuDd+L06FQ3h5mERDYq/G4PdDdTXs2AHnnQc/+QlMmeK2CgyHIXHkWfD883DEEXDBBfDRj0Ig4LaUCIfdajMiIqNQt1pEGGMOt9Y+PwDj6RNqESEiIjKwFt62kLUb1/b4vPIF5UQWRvb6/tbCzJlu0g/c6aHBIMRiuL0Df/AD+Na3Wk9wHDcgjEbdnoKVle4J0agbXYqIDHOdFYbpbhB4CHALsI+19lBjzGHAR6213+3bofYNBYEiIiKDp7cBIXQvKIxEoKKi59duwkMB7dcvth+EKoaKyPDXF0HgWuBK4GfW2ne1bHsp3fZhqFEQKCIiMvgSDQnKVpRRWV9J0BfknT3vsOntTQCZHoTxS+J9dr/sBvNpV1wBP5xwDQXf6f4axrwUGIrIMLM3LSLSxltrn2yzrXnvhiUiIiIjmX+an9jiGKnyFLHFMdZcsIa5M+YCYLFMHzedN7a/0Wf3i0bdGZ7pqaDnnQc33AAn/utq3njdQm2tW0LUGDjoIHf+aDZHrSREZHTobhC4xRhTQks/v5aG76/126hERERkxPFP81N1aRWpq1P8ouwXPPXqUxzw4wMwFabPCsXEYm4mMB6HO+5wH48/DocfDqs3ZB2wcaO7DvATn2i9wKxZaiUhIqNCd4PAS4CfAQFjzH+BLwEX99egREREZOQyxvC5Iz7H/pP3pznlTiyqrK/k1N+d2uf3eu973YTfG2/A8cfDkiVunRgAJk+Gu++Gv/zFPeiVV+C226Cxsc/HISIylHRrTWDmYGMm4AaOu4CzrbV39tfA9obWBIqIiAwdkTURKtb2opILe189NN86wYUL4c4728wGbWiAL33JTR2OGQN79rhTQ1UtVESGqV4XhjHGTMbNAu4P/BFY1fL6K8Dz1trT+n64e09BoIiIyNAXXh6maksVKZvCYADYf/L+3H3m3Rx34HHdvk5vq4Vmy9R9Oegg2LSpdYffD2PHqoWEiAw7e1MY5jfAXOBF4CLgIeAs4PShGgCKiIjI8BBdFCVQFMBgCPqCPPDJB/B6vCy4bQE3PHYD3Z2tFInkr+VirZvMc1p+23Ec8HrBZLWJmEE9TRQQqTBuwZjsABDchvLxuHuxeBxKStzjVChGRIaxrjKBL1pr57U89wBbgIOste8M0Ph6RZlAERGR4emt3W/xmT9+hvur7uf0wOn8+rRfM3Xs1B5fp7fZwRJq+CunMIfanp+sNhIiMoTsTSawKf3EWpsE1g/1AFBERESGryljp3DfJ+7jxyf/mJUvr+SInx3B068+3ePrdJQdbJsZ3H//1tcA650STvf+zc32hUKtmb80n8+dFpreX1urNhIiMux0FQQebox5u+XxDnBY+rkx5u2BGKCIiIiMDomGBOHlYTzXePj5Mz9nxcdX0JRq4rhfHcdPn/ppt6eH5hOJuHFb/P/bu/PwuMq6/+Pve5Ju0AKFVvYtBZom8ACyCjwUBRHFiLJpQRZBQcuOgIBApuiPVR5lK4IsCmLZUYIgKFpAlH0RkqZQRoQCQlNaoKWlbeb8/jjZO0kmy2Qmmffrus6VWc45cyfAlXy47/v7rWstEpNOw9tvty8ak04H6pZtRojShLpayl5/mNpoEitIsGzYKjBvXrw/MIriijNVVX37piUpD7oMgVEUlURRtFrTMSaKotI2j1cbqEFKkqShr2pGFfUN9URE1DfUUz2zmheOfYE9N92TH/zxBxx6z6F8/GnvFiR1t2+w7WRfaaKRTUgRSDOKJYxiCaWkGb78k/Y3bW5IGELmw5lBSQWqRy0iBgv3BEqSVLj60jJirVFrMfPImWz5mS37bTypVDyhN2tW3C9+7lxYEbcvJIR49WdtbdPJlZWtM4EAJSXwk5/ALbfEM4NWEJVUIHrdImKwMgRKkjT4tG0ZkQgJShOlrEivaHm+4WobsnTFUj769COu2fcajtjmiF5/Vm8LxwxjGfdRxV48wn/YiFHrrsl677bZs5hIQHl5m9QoSflhCJQkSQWrLzODZ+92Nv9vz//Xr+NZucF8RCAiIkGCRsqpp5ZezkRaQVTSADEESpKkQaermcFmW31mK+46+C62WGuLPn1WXxrOv0Il5cyihIgI4rb3IcRLRp0ZlJQnfWkRIUmSNKCSM5OEaYG6eXUtgS8dpVnWuKxdAAR4+f2XmXjVRMK0QJgWSM5M9u4zk50XjmnuCNFs7bVhxbnTmuYHA5XUUUL8P9VbTmv+n+wWj5FUgAyBkiSpoCT3SBJVRy1H9eRqouqIivEVJEL8p0siJKgYX8GbJ7/J5zb4HADH73A8Z+12Vr+P5/774yAIMHo0vPcefPYP1Uwoi0iEiMqKiNTrXZQbHTECNtyw9bVEIj7H3oKS8sTloJIkaVBILUhRNaOKWfNmMWn8JK7+ytUc98Bx1M2rY81Ra/LBkg/YYb0duOOgO9hkjU1yMoYogttvh29/Gxob49dCgEnj3qd23tr98yHuG5TUD1wOKkmSBr2ysWXUTq0lXZ2mdmotxz1wHPUN9QAsXLqQDVbbgFfnv8q2127L/a/e3+fPa24w3/ZIJGDKlNYACHEwrJv3mabFoe2PZHUEy5bBuuu2v/mmm8Lrr8czgiUl8dfXX3d2UNKAcCZQkiQVtL5UDwWonlxNco9k/w2ITBVEW2WsBZNKwb77wuzZcZpMp2G11eDjjy0gIyknnAmUJEmDVsc9gs1Hpj2CUXXE2budzVHbHAXAXmV7cdwOx/XPOJKtM4J1dZkDIHRSC2ZCGaF+FiFKc1H69PjEjz6ygIykvDAESpKkQalmSg3l48oJBMrHlXP1V66mcnolF/z9Ap58+0ku3PNCHv/P43z2us/y5Nwn+/x5XVUQ7VgLZpVV4MEHW19vu9rzzOii+MEWHdparLYalJXFs4JgARlJOeNyUEmSNCR07CtYPq6cW75xCwfecSBzP5rLz7/0c6buMJXQNq31k1QKqqriybx11oElS+DDD1vfTySgfK1+LB4DFpCR1CWbxUuSpCGjr3sEz9rtLC7Y84K+jyPZ+wbzHVVXQ/LED2DiRGhoaH1j443hr3+NE+bs2fH7NTXxjKEkdcEQKEmShrxMM4G1U+NCK+kozYWPX8i5fzuXivEV3H3w3UwcNzG346mEWbNat/2NGBHPEr75Ztx3MGOWS6XgK1+JA19JSVyGdJVV4qlFC8hI6gELw0iSpCGv4x7Bmik1Le8lQoIf7/5jHvr2Q7y3+D22/9X23FV3V27HU9PaZH799eM895//xFlu1qx4cm+lPhQTJsQBEFr7UHzyiQVkJPUrZwIlSdKQ1bHBfM2UGoYlhnHQnQfx1NtPcerOp3LRXhcxrGRYv3xevy8RTRInydmzW4NgSQlMnQoPPQSvvdbFtKKkYuZyUEmSVJQ6WyJ67l/PZeHShVz1zFXsttFu3H7g7aw3Zr3cjiVDb8ExY+CMM+CUU2DVVTu5sLnqTH09bLghbLYZPPJI6/shxEHQJaKS2jAESpKkIauvhWIAVh22Kn885I9M3mRyP41qZc1ZbtasOLNdfjlceik8/HD8/jrrxDVgJt2e7L/pxGZWEpWKjiFQkiQVpY4zgaWJUlakV7SbGbzjwDvY/479ef2D17lwzws5bZfTctJGIuP4OhSPKS2F6dPhyCNhWFcrVNtOK4YQX7h8eev7zg5KRc/CMJIkqagkZyYJ0wJ18+pIR/H6y3SUZlnjsnbP6+bVseU1W/Lq/FdpjBo54y9nkDg/wZl/ObN/xpHsvIZLCHGNl7b/P37FCjjmGBg+vJO6LyHZemHzutIoah8Am1/rqoCMRWSkomYIlCRJQ05yjyRRddRyVE+uJqqOqBhfQSLEf/4kQoKK8RUt55y3+3n8397/R0ko4e5Zd/Ov9/7V93Ek4zzW2VFREXd9gPjrpElw//2wzTbxa5tsEs8MNneISEZd3LDtzZpvCHFvihDiHoOvv956viFQKlqGQEmSNOQl90gCmdtIpBakqJxeyfmPnc/1L1zPrfvfyuJli9n5+p255aVbcjqumpq47V8I8df774d994Xnn4/fW2eduBDoppvCz34GH3+cxc1KSuJA+OKLsO668OmnceibPRu22w6efrr99KOkouOeQEmSVNQyVRB95PBH+NZd3+LR/zzK97f7Pr/Y5xeMKB0xIONpW0CmvBzOPht+/eu4IOjYsXDSSXDCCbDmFcn+LyADFpGRhggLw0iSpKLWHxVET9rpJH6xzy/6Z0BdaFvzJZGIg2BNDey5J7zxRnzOqqvCoYfCzJndtArseLPNN49T5CmnxDOEzTbZJE6fA1QQR1LuGQIlSZI60XEmsGxsGcNLhrc0mD9+x+M58y9nUpooZcYBM9h7wt59+rz+bCjfrF0x0Fx8ADhDKA0yVgeVJEnqRMd9ggD1DfVERNQ31HPV01fx7PeeZb0x6/Gl336JCx6/oKXCaG9kKhZTXd2376FdMdBpSQLRSkeyupMiMiHEewd33rn9TZuTpYVkpCHHmUBJklSU+mOJaPXk6paiM/2lY1P5mpr4ecclotD6GsT1YBob4xWgJ58cLxcdNarpprmaHQRnCKUC5XJQSZKkLGUqFAO0vBYIRERssdYW3PvNe6kYX5HzMTUHw9mz404PNTXx623D4l13xa+fd1683a+kBI49Fn78Y1hvvSw/qGMT+lVWgU8+aV9N1Eb00qBgCJQkSepGoc4M9kRlZRwK2/55N2wYfPOb8ezgdtt1uMAZQmnIck+gJElSNzprMJ+pyfxma27G5mtu3nLtdz/7XZafuzznATCZbNr318lRV7dyC8Dly+G3v4Xtt89wzbRk+72CXTWib1s5tPnxDjvAL38JH37Y/vzXX4c774zXrlZWxlOZkgqGIVCSJCmDtoGuuXgM0PL19QWvt7x//fPXs/cte/P+4vdzO6Zk91mtud5LIhE/j6K4b/zaa7feZ+xYOOOMOKtlNVlXUxMvAW1uRP/MM/CLX8Th7/vfh9VXhzXWgNtvjz+weRNjFMVfq6r6+0chqQ9cDipJkpRBfywPhYFdIpqpqExZ2cpb/UaPjrf6pdPxtr8lS+K9hvffH9+n4/7DsjJyu3S0mUtIpX7jnkBJkqR+lKl4zC3fuIX9b9+fdxe9y1Vfvorvbfe9ARlLLrNZgkbKqaeWLUmxKVXUMJuJTGQ2NVRRxr97d2PDnpRz7gmUJEnqRx2Xh9ZMqeGz636WgysPZo9N9uCY+4/hu/d9l6UrluZ8LLnoO9gsTQl1VBKImECKOipppJQ6KplAiloqaGz6czId2qw/ffhhWGut1hutvTacfTbMmWPPQakAGAIlSZJ6qGxsGbVTa6meXE3t1LhVQuX0Si79x6W89eFb/GD7H3DDCzfwvzf9L29++OaAj69jMHz99dbaLhUV8fNM9V56qooa6ilnBSXMisopq6uJC87s/UXC/AaGsYzbDr4nrkpz0UWw2WZxU/oNNog/uGPRmFQqfq201IIyUg65HFSSJKmPMi0PveALF3DYvYcxonQEtx1wG3uW7ZnvYa4k0x5CgH33jfcQjhwJS5fGTedHjGgtAtrcsL5HrQLnzoXf/CZeu7p8eevrm27aGvbabl7s1YdIauaeQEmSpD4ajIVi+sPzz8ddIH7727iADMC668YrPrfcMn6eOulyqq7Yq3/2C/aW+wyldgyBkiRJOZRpJrB5meiiZYs46g9HcWfdnRww6QBu2u8mxowYk+cR99yHH8Ktt8I118Arr8Bqq8Hhh8Oxx8bN6LOewOtYqnTttePjpZfi581/m4YQT086Eyj1ioVhJEmScihToZhmP/vHz7j9wNv52Rd/xr3197Lj9TtS31Cfr6FmrWNj+jXWgOOOiwMgwEcfwVVXwVZbxU3q0+n49XQ6ft5ZQ/uyuhpq0+WkQ0kc8p54Al54AZ59Fg45pLXR4YgRcMQR8XpUSf3KEChJktRHHQvFlI0tI7UgReX0SqY9Oo0tr9mSb0z6Bn8+7M80fNLAjr/akXtn3ZvvYXepu8b0UQTz5sGll8KwYe2v3WyzzpvXp6IyKqNaEukV8SxfWVmcDrfbLl5z+tFHcO21MGEC/OhHsOGGcM458I9/xLOIiUTmojHNRWU6e19SC5eDSpIk5UBnS0Tf+vAtDrzzQJ5++2nO3PVMfvqFn1KSKMn3cLM2YD3jqyP429/giivgvvtal4lC5jWnFpWR2nFPoCRJUj/rr0IxAKfvcjqXfPGSfrlXvr39NtxwA1x3Xfy42Upb/AYiTWbLojIaggyBkiRJA6zjTGDZ2DKGlwynbl4dFeMrqJlSw9/+/TeOe+A41h69NvccfA/brbddvofdb1asiMPgD38IixfHr33ta3DWWbDTTr3oT1hZGfeyaPu3a0UFnH56vJdw222dCZTasDCMJEnSAOtYLAZoKQhT31BP1Ywqjv7s0fz9qL+zcOlCdr1xV2564aa8jbe/lZbGlUMXLYJXX4WTToKZM+Fzn4Mdd4TLLotnBjNt4cu4va+mJr6gpKmgzM9/Hm9G/M534l6DVVWw+ebx++XlrU0PJa3EmUBJkqR+0p9LRAdbP8FsfPwx3HxzXFW0vk2B1I5LRbPe3hdF8Oc/wyWXwCOPxH0rTjoJTj4Z1lxzIL4lqWC5HFSSJCnPMhWKAVpeCwTWHLUm85fMZ+cNduaug+5i/dXWz/Oos5fvLX6f5TnO4kIO5G4YMyYOg6ecEofBVCqeKZw9GyZOjGcJy8ryN1hpABgCJUmS8sCZwc61ne2D1j7x++wT9xl86634eceZwG5nCV9+Gc4/H+66qzUM3nknvPaa+wVVVAyBkiRJBSJMC0TV8d9fnbWRqJtXx+437c7CpQv52d4/46SdTiL0uJJKYWuenOtY66W/bMnLnMf5HMRdRECff3pWENUgY2EYSZKkAlE9ubrlccfiMTVTakgtSHHQnQcxf8l8Vhm2Cqc8dAqH3nMoi5ctzteQc6KsLJ6MS6fbN6BfuhRuugm23DI+b+ON4cor4wqjnTWgz9TI/uVoKw6K7oR//YswZkzrB4cQzwRmuqirwwCoIcQQKEmSNIDaLuksG1tG7dRaqidXUzu1lrKxZVTNqGqpIrpo2SLGrzKe2165jc/d8DlOfPDEPI164IwYAUceCf/6V7x1b/314YQTYJNN4IIL4NZb4wzXnOXaFgHNWFV0q63gxRdb9wBGESxcGCfNxsYB//6kQuByUEmSpAHUn/sEYejtFczk8cfhwgvhwQfjbX4/+EFc82Wdddqfl1VV0ccei3sLPv10PN14ySWwxRZxE8NZs+IypRaO0RDgnkBJkqRBorMqorPmzSIi/rtt/Crj+e9p/yURhu6irtxWG404kLu4kLPYjNdZxCqMYgklZKhEIw1ShkBJkqQC5uxgz8yZAz/9KdxyS7x8dOpUOOMM+PyW71M/by3SlJCgkXLqqWXLgRmUhWNUYAyBkiRJg0xnVUQDgbVHr03DJw1svPrG3PvNe9lq7a3yPNr8eO01+MlP4n2CI0fCIYfEqz1fe62HqzrLy+Megs3WXjvuUTFsWM7GLuWa1UElSZIGmc6qiE4aP4knjnqCmUfMZPHyxex8w87c9spt+RpmXm2+Odx8c9xXcP/94cYbYe7ceMvf44+vHAAzFo4BeOCBuMwowOjR8N57sO22caKUhiBDoCRJUgHqrororhvtyvPHPM+262zLlLuncOpDp7K8cXn+BpxHEyfGS0Nra+HrX4dLL40D4MUXw5IlredVVcWFY6Io/lpV1fRGc7+KKIKPPoLf/x4WLYLJk+Gww+C//+0iQUqDj8tBJUmSBrFljcs47eHTuPLpK5m88WRuP/B21h69dr6HlXO5LRwDo/iEs7mA07mUpYzkE1bhM7xn8RgNGi4HlSRJGqKGlwznii9fwc1fv5mn3n6K7a7bjifnPpnvYeVcMtl9f/eZM2GnneLzJ02CDVdbSCDuDZigkQpqiQgZj09YlXP4f4xgGavzEevy3zgAQtyDoq4ublbYk8PCMSoQhkBJkqQh4LCtD+MfR/2DYSXD2P2m3bn22WsZiiu+emLyZPjnP+Huu+O+8G99tAYjR5UAUF5RQs3rld0nySiKQ98GG7S/eXl5dte2PQyBKhCGQEmSpEEsOTPZ8njbdbfluWOe4wubfoHv//H7fPe+77J0xdL8Da4AhBAXjXnlFbj2Wlh99fj1//mflYt/drrtLwR49NF482GzxkZ47rkB+R6k/mYIlCRJGsQ69hdcc9Sa/PGQP3LO/57DjS/eyP/e9L+8+eGbeRpd4Rg2DI45Ju4xeN55ce2XiRPh/PPhk0/iczotHANx8ZjmNx94ABYvjteaVlfDsmXxORaP0SBhYRhJkqRBKLUgRdWMKurm1VExvoKaKTWUjW3fE+EP9X/gsHsPY0TpCG474DYef/PxId1Evlmui8YArMECLuckDucWXmAbjuTX/I5DmBTqSURpi8co7ywMI0mSNMRUzaiivqEegPqGeqpmVK10zn7l+/HM957hM6t+hr1/uzfTHp1GOkpndf+2y0wHm2yKxrQtHrP11vF1qwxb1m3hmNcpo4JaPmR1LuZH/JfPsC0v8hLbUEldHACh98VjLCCjAeBMoCRJUoFJzkyutMyzr0aUjODTxk8ZPXw0jx35GNuuu23LZ2WaHQzTAlH10Ps7MZPGRrj+ejjrLFiwIH5t4sR41WfHhvOVlfGq0HTbyb5HG+D44+H221tPdCZQedbVTKAhUJIkaRCqnF5JfUM96ShNIiQoH1dOzZSajEtEK6dXMmveLKKmFgfDEsO48+A7OfuRs1c6N5tlpoNdrpaLTuUqLudkSmjkbdZjdx7n36z8s6uudrJPuWcIlCRJGmLahrX+FghERC3hsnZqcc1m/etfcOyx8OST8PnPwzXXtBYGzTgT2PbH8847MGUKPPYYHH00XHkljBqVl+9Dxc0QKEmSNMjlYolof6meXD3kCs6k0/CrX8GZZ8bVQ888M14u+s47cdXQujqoqICampWXjLJiRVyC9MILYaut4M47W1NkKhXfYNasuIN9xhtIfWcIlCRJGqLa7t3LtES0dmptxiWeX/3dV6lvqG9ZIjph7ATmnDin03sUq/feg1NPhd/9DjbfPJ4V3HPPeDlnt0s6H3wQDjsMPv0UrrsuniHsdipR6h9WB5UkSRqiqidXtzyumVJD+bhygJY9ggBlY8taglzt1FrKxpZx/yH3M2n8JABKQgnvLnqX2165rdN7FKu114Zbb4WHH46rie61F3w7/Jap0z7TfZXPr3wF5s+HRYvgkEPi1+rq4gAIfasgahVR9YEhUJIkaRBruwwzU9hrq21gbHvuGye/wTbrbMOUu6fw83/+nOePeZ7qydUZ71GsvvhFePllOPdcuGPYt5m05vvc+tuIKJ1FH4ply+D00+MbjRwZBziIZwIrKrLvZ9HxMASqlwyBkiRJRSLTvr3qydVssNoGzDxiJqfufCpXPXMVu/96d76zzXcGfoAFbuRIOP98eOkl2GIL+Pa3Yb/94n2CmbRktGHD4JJL4L77YMSI1hBYXh7vCZQGmCFQkiRpiGk749ed5mA4rGQYl33pMu4++G7qG+rZ9tptuW/2fSufP4ibyPeXSZPg73+H//s/+Mtf4sm8m26KJ+cgrv1SWRm3oaisjJ8DcUGYl16C7Zu2aR14IMnfbJqX70HFzRAoSZI0xPSlUuf+k/bnuWOeY5M1NmG/2/bjlD+dwrLGZS3vd1ahtNjCYUkJnHJK3E5i663hqKPiLYBvvhlnvfr6+Lz6+vh5i403hscf5+MDjoTzz6f8/Cl8dtKS1qAoDQBDoCRJktrZbM3N+OfR/+SEHU/gF0/9gl1v3JW//ftvVE6vBOIqpKkF7VNLobav6I1kMvvaLJtvHrcEBPjTn+KM123tlxHDWe3uG/kRF/Etbmd6/efZZcJ/rfuiAWMIlCRJ0kpGlI7gii9fwb3fvJc5H8xhr1v2Yta8WQDUN9RTNSOe3kotSHUZDgejZLL7mizV2a+47UTgEn7E/tzN//AvnmIntuJf7c6YNs0iocoN+wRKkiQVsVw3oR+KjeS7EkVxj/hzz41nAdddFx5/HCZMaH9e23aB24fnuL/ka6w98iO47TbYd9/8DF5Dis3iJUmS1COGw755803YYw/497/jTPerX8WBsFkqFe8VrKuLC8s88Ku32fiEr8GLL8YVZ048sbWKqNQLhkBJkiT1SWpBiqoZVdTNq6M0UUqCBBfudSHXP389s+fPJh2lSYQE5ePKW/oPFrt0Gq6+Gs44A1ZZBX75SzjooPbnJJNtlnMuXgyHHQb33gvf/z5ccUXcXkLqha5CoHsCJUmS1K22zeX/+8P/8pUtvsIPH/4h41YZx4Sx8VrH8nHl1Eyx712zRAJOOAFeeCFeDnrwwXDoobBgQes57fbzrboq3HUX/OhHcWLcd19YuHCAR61iYAiUJElS1qonV7PWKmtxz8H38Mt9f8mz7zzLgqVxqqmdWkvZ2LKVrunv9hGDrR1FeTn84x9xo/k77oCttoKHH+7k5EQCLroobjw4cybssgvMnTuQw1URMARKkiQpa837+EIIHLv9sTx3zHNssNoGAHzvvu/x0acfrXRNf+8tHIztKEpL42IxTz4Jq68OX/oSHHdcvAI0oyOPjJPi66/DbrvBnDkDOVwNcYZASZIk9dqk8ZN48ugn2XXDXbnxxRv5n2v+h5lvzMz3sArWdtvBc8/BqafCNdfANtvAP/+Z4cRUKk6Jy5bFM4Gf+1zcmb4r9opQlgyBkiRJ6pMRpSP4+1F/5/HvPM6wkmF8/jef5+Q/nUzdvLpe9xDMtORzqPQkHDkSLrsM/vpXWL48nuj78Y/jvNeiqiruIQFxhZmFC2Hy5HgqsaNUKu45MW1a/DU1OH8uGjhWB5UkSVK/WbxsMWf+5UyueuYqhpcMZ3njciKiHlcODdMCUXX7v1Mrp1dS31Bf+JVIk8k4kGXhI8ZwCj/nRo5ma17kFg5jK17p/zFVVxf1TGFyZnJItyTJxBYRkiRJ6he57h+YSwXVm7BDULyPKr7Hr1jIGvyEc/khl1FCOj9jG4KBMdP/VBjqDIGSJEkacOVXlTN7/uyW5xuvvjFvnPxGy/NMf5i37UdYMb6Cmik1LRVHB81MYC/Nmxe3B7znHth1u6X85sOvM2HOQ3E3+ZoaKCuDDz6IW0c88wzceCMcfni8BLS+Pl42mkjE5Uhrh87PpT8YAttzT6AkSZJy4oFDH6BifAUAwxLD+M+H/+G7932XBUsWdHpN1Ywq6hvivXD1DfVUzahqea9mSg3l48qBodmTcPz4uE3gLbfAK3NGsvW7f+I6vkf0Sm0cAAHWXBP+/GfYYw844gi48so4IJbHPxfKy+Pn6heDrR1JtpwJlCRJUk6FaWFAPqeglnv20VtvwXe+A488Al/+MkyaFBeTabF0KUyZAr//fdyA8Jxz4lnAIfi3fX/o7UzgYJ5B7GomsHSgByNJkqTi0/yH9AvvvsAuN+7C0hVLAQgEJo2f1LKsM5sln4P5D/Nsbbhh3CZw+nQ44wx48MG4ncS3vw0hEJcYvfNOOPpoOO+8uHroeedl/wHJ5JDb96fsORMoSZKkPhuogjHVk6vjzxsiM37dSaVg773jnvEAu+8eLxfdaKOmE9JpOPnkeFnoD34AV1/dlBK7EUJRzRo6E9ieM4GSJEnqs+Qe2Zfgbzvb12yj1Tfisr0v44BJBxBC6NMf30OpHUBVFfz73/HjEODxx+M6MBddFGe+RCIBl18ezwxeemm8JPTKKzsPgqlUfFOIb9RccEYFJ5f/HlsYRpIkSQOqbYGXivEVzDhgBmuMXIOD7jyIPW/ek5ffe7lP9y+EFhbJmUmSyTiL9eWoq4sn+yCeuIsiWLQIjj8eSkqazksEwqUXcxmnwtVXc3niZEKIMt7v/Z3bNKGvr28NhENUakGKyumVQPw/H1ILUnkeUfZy+e+xIVCSJEkDqmxsWcs+v9qptXxry2/x3DHPMf0r03npvZfY9tptAXh/8fv5HGafTHt0Gslka3Dr7VEx/n0SNAKQoJEKakkT+DVHMJYPGMFSLuAsljGcH/J/AJzEFUQkiAgrHZ+Z1yZVptNxyuxrUi3gvYVdVZstZu4JlCRJUl5kWvI5/5P5nPe387jm2WtYZdgqnLzzyZy2y2msMXKNPt13oPXXGJpXb9bVtW8XCPDf/8IJJ8RtJbbZBm64AT67bcRTO5/ETk9fCaedBpdc0n5paAH0FByo/aP50J8Vavv675DN4iVJklRwuvojd3bDbKpnVnN77e2sMXINztjlDE7c6URWHb5qn+47UPp7DF3VcbnnHjjuuLjZ/Omnw0UXRUTHnRAXifnRj+DCC1uDYFepcgjKptpsV/L571IuQ6DLQSVJkpQXzZU+M5k4biK3HXgbLxz7ArtttBtn//VsJlwxgSufupJPV3ya8ZpC2P+VjzHsv3+c6Y44Ii4YA4Hbdr2S9LE/gIsvjnsINifIsrKWmb/k1QcN6QAI7feflo8rp2ZKTZ5HVBgMgZIkScqLbJbNbbPONtRMqeGJo56gfFw5J/7pRLa4agsuf/JyPv7043bnZrv/Kzmz+8/trVztQavuPC8DMHYs/PjHsPHG8fMphwS2+cfVzNz3UrjggriHYNupxOrqIbsks62O+0/Lxg7t0Jstl4NKkiSp4ORj31jzfq5C2LPWm71lbbf7hRBXD12xAvbd8CUufusQKqsPalfEpUfLDQd5c/nB2CfQPYE9ZAiUJEka2p6a+xSX/fMy7p51N4mQYMqWU3jirSd4Y+Eb3e7/yuUf9j3dg5ZMwrSc5824NuhR3MhoPuZyTmn68ADJrn8O1dVN2W+QN5fv6T/z1IIUVTOqqJtXR8X4Cmqm1Az4LKJ7AiVJkqQ2dtpgJ+446A7mnDCHqdtP5Z5Z95BakGJU6SgAJq41MS/7v3q6By3bNhLVf0t2306iIi74CfHXior49YaGwMknwc2JI7mOYzlnj7/z4Yfxed3dcxBP/vXJUG8t4UygJEmSBr0FSxZw3XPXccXTV/DOx+8wfpXxHLLVIRyx9RFss842hDZtErqbYUnOTPa5zH+/VwfN4n7dFf5MvdbIOV98ihn/2YXRwxtYtFc14z77d5469d7uZ7nyNBNYCEtzB0rHJcAuB+0hQ6AkSVJxWt64nOE/Hc4Bkw6g5tUaljUuY6vPbMURWx/Bof9zKOuMXqfbP677I8DlIwS2nNtVXmts5Bdb7ckfZlUzk8/DGv9m/G41PH/NiWywQYbzh0hLiZ7+8+hra4n+kMsQWNrru0qSJEkFZljJMADuOvguPljyAbe9chs3v3Qzp/35NH70lx8xaljTctGrJvLgoQ/mvVpkNrOOzdVMw7TQ5XmtF0DoYvKs5EC45e4vcFrdV/hB4lTeuv9ENvxjI0x4CD57A2xRA6XLAXjlaihvgBKgcVYd9TtPYMvjVr5nfzZJLwQ1U2pa9gQOxdYSzgRKkiRpSMk0g1LfUM9uN+7G/CXzW14bWTqSM3Y5g69s/hV2WH8HEiHR6fX9MYacn9dUZSYQEdF5YKycCq+tCbfeAwfVwenbbcrI577DrzmSuWzIOOZxGLdwNDdQSV23Y+uVloozA8PqoO05EyhJkqRBJZt9YtnMmi1dsZTzHzuf8x87P6vrezLbVT25m8Z+uZBMktwDqmfSZdXPmgUpvvSbKg75Rj2Jpaty6XP/pqF6LMlzN+Thh+HGG8dz1R9O5efLT2XHkf/i6E+n863od6yWWAzl5S3N5jV4ORMoSZKkIaWzJZaZ9nk9duRjPPT6Qzzw2gM8OOdBPljyAQCbrLEJO62/U3xssBPbrrNty1LS/tTdbE9PWxVkO3tUWQl1sxoZFjVyJwezH3+Aq6+GqVMBmDcPfvtbuOGaZdS+NpxVWMxBq/+ZQy7fiZ2+vi6rr97z7zWfnAns8J4hUJIkScWgu0D12vzX2Pu3e/PGwjcYM3wMo4eP5t1F7wJQmihl67W3Zqf1d2KH9Xdg0rhJbLHWFowdNbZPY+ruD/22wZV0AhrKYXoXM3FZ9P7raBjLuIsD+Ro1HMsvuY5je3R9WwO8yjNrhsAO7xkCJUmSVEw6++M600zhXw77C0+9/RRPzX2Kp95+imfeeYZFyxa1XDNulXFssdYWTFxrIlustUXL49+9/Dsu+PsFA/ltdavtctbmmUCiEhIJ2GqLT3lxwgHwxz/CddfB977X/uIQWLwo4u9/h2eeiY+nnoL33ovfHjYMtt4adtwRdtgh/jpxIpSUDOz32Jnetv3IWwhMJglhmiGwJwyBkiRJxatQesuNHTmWdcesy3pj1mPd0evGx5j463pj1mPdMeuy+ZWbs+LcFZQkMqelnrYqyDa0pFIwYedamFfZ2vlh/U/hG9+ABx+EG26Ao45qc+OV+05EEcydC08/HYfCp5+GZ5+Fjz+O3x8zBjbZBNZaKz7Gjev68Zgx8ce0PbKRTsOKFdkd6XTrt9N8/0yf2fx48ysnMOfE1zt9v7PX2j5PJGD8+Oy+l7Y/75DEENgThkBJkiR1piczgT3pDffxpx/z6vxXeXX+q6QWpHh30bu88/E7vLvoXd79+F3eXfQuyxqXZbx29PDRrDVqLdYYuQZjR42Nv44cSyBQ82oN8z6Zx7qj1+WMXc9gwtgJjB4+OuMx/KfDs+8nOC1eOtouDixdCl//Ojz8MNx0ExxxRNPJ2TWLT6dh9uzWYDh3LsyfHx8NDfDBB9DYmNXw2o+1k3C4YkXP7zXQRo9uDcZZy3EItDqoJEmSRN97w40ZMYbt1tuO7dbbLuP7URTxwZIPeHfRu+z7u31568O3iIj/yB9eMpzdN96dhUsXsmDpAuZ8MIcFSxawcOlCFi9fDMC7i97llIdO6XYca12yVqchcfSw+Ouqw1eNT95hOje/1OGca89n9PGLGH3ckYxmOaMOP5pQnV2100QCJk2Kj+b82FY6DR991BoK2wbERYvinNn2iH9umV8DKC3t2ZFIZHfv5uPwew/jN1+/pUfXdHxeWoCJy5lASZIkFZXulkx2936hLDcFOGDSAeyz2T78Z+F/uPa5a5n3yTzGjhzL5zf9PIHAomWLVjoWL1/cbl9jdwKh81DZ5lh12KpZndccQksTfUtHyWTui9DkbU+gM4GSJElS4Uju0bsiI231916/yumVzF8yH4APP/2Q+ob6Lu+XjtKc89dzuHC/k5jz5spBcdGyRSxa9AGLrruKRW+nWPTlnVi0+SYsWt76/rxP5vHGwjdann+87GNWpLNfnzmydGSvA+Wo0lFM+/UovvK9kYwqHcWoYaNavo4sHcnI0pEkQiLrsRQbQ6AkSZLUA72tNNlW26Wn6ShN3by6bhvcd/d+s2zvB8DpF7LZlSu/3FJJdNtj4IAD4PSH4Mc/hp9c12XFlmWNyzIHym6O5tnJRcsW8f7i99u998nyTzJ/2Hdgp+s7/9ZGlIxoCYcjS0cyvGQ4w0qGxV8Tw7p8PjwRPwY4+U8nkwgJEiFBSShpedzd8efUn/nyZl/O6tySRIf7VnT/j64vXA4qSZKkotJdiOvrctGeyPZePekn2JOiNlktqVy+PG4if/318O1vx5VDhw/v9t5Z378bjelGPln+SbtguGTFEnadvISaPy1h6YqlLFm+hCUrlrR87fja0salLG9czrLGZSxPN33t8DzTaw2fNLD6iNVpjBpJR+mMRy5l6meZLauDSpIkSVkajCEwtSDVMrPYl+DQqSiCCy6Ac86BPfaAe+6BsWO7H3d2RUV7JZf3bvmMLP75RFGUMRyOvnA0C3+0sNPwmI7SGcPlV3/3VVLzXydK0KsqtS1jd0+gJEmSNHSVjS2jdmotYVroVWDoVgjxctBNNoHvfAd23TXuJ7jxxiufOxAVW/pZV8V+sl2Gm8kaF6/RuwubtjP2aGlvDxgCJUmSpAJXPTm7Fg05d+ihsN56cVP5nXeG+++H7Tq0xJg2bfCFwE6K/fRl1re311ZOr6T+vTrSfZ0JTHYeHC2ZI0mSJBW4vhai6ZcxzGwaw+c/D088Ee8LnDwZHnggr+Maamqm1FDeED/uTb/KbBgCJUmSpDYKZtatwLRbLllZCU8+CRMnQlUV/PKX+RtYdwbZrGTZ2DJqp8ePa6fW9u/eziaGQEmSJKmNQph1GxTWXRcefRS+/GX4wQ/gzDMhndtqmb0yLfNev2JmCJQkSZKUlZUm1UaPht//Hr7/fbj44njPoAqeIVCSJEnKk8G29DTjpFppKUyfHofA226LX3v66QEdVyErxNWohkBJkiQpT4bM0tMQ4Iwz4K674uc77QSHH876zM3vuPqoP0J6Ia5GNQRKkiRJWUgtSFE5vRKIy/inFqTyPKICdMAB8dezzoI77mA2E+MU9Mkn+R1XLw2ZkN6BIVCSJEnKQtWMKuob6gGob6inakZVnke0skJZXpr67gXsveEs/si+kEyyYsJEuPXWfikck0rFxUkh/prKVRYvxHWc/SREUe+aHxay7bffPnr22WfzPQxJkiQVoOTMZPt2BwOsenL1oJxhCtMCJCO6jQ8hUFkRUV8fZ77dw+NMH3EylUufj5eJ/uIXcaP5XqqshFmz0kRRgkQCysuhtqte6iHQ/aD78brmy5uaxffqNiEQkvS6UX18i/BcFEXbZ3qvtNd3lSRJkgah5B7JXoWwyumV1DfUk47SJEKC8nHl1E7tKn0UtmSyh/vVkvGXEDK/vSkpaqhiIiXMrltBuilqPBb9L1stfYbDuZkLnjqb9T73OW7lEM7kIuayYS9HHy9oTKehrq7zMcUi6PL9zKqpbv6WhxyXg0qSJElZqJlSQ/m4cgDKx5VTM6UmzyPqm2QynqHK9mjW2fupiioqE/WU0shEZpOgEYBEAiZVJPh1dCTrffwqnHMOh468h7dGTSQ6/QyiR/5K9PEiqquzG0dFBRBa711R0c01hB59n81HkgKs6NJPDIGSJElSFsrGlrXM/NVOraVsbFmeR9RLyWQ8ddbTo1ln79fVtez5q6GKcuL9k+XpWmrqyuJzxoyBn/4Uli6FJUvg0kthzz1hzBi+Pm2bzu/dZn9eTQ0wrune5U3P1SOGQEmSJKmY9HQKMNupwIqKeGoOKEv8h9qKgwGojSopi1KZr/ngA3jwQTjvPN7nM3FIbDZ+PHzta3DhhfCFL7RUGC0rA47bMr53bdPzQlPgRWXcEyhJkiSp72pqoKoqnhFsnqKb0M01Y8fCPvvAPvvwpfMhWtAYX//Pf8I//hF/ve+++NzS0jjxDRvGS/MgzTbw2UQcPLs6IA6RsHKYbavt8+aZz9NOg402go03bv06dmx3mxDjzZbJbr73PDIESpIkSepSKgVc/QoQV+esqckwA1dWFk/NhdBNuc4ulJTAVlvFxzHHxK81NMCTT8aB8LXXIJ1mTmMtifc3YZsN0vES1M6O5cvje6xY0foZHZe2tn3c/Ly5lcX06fGy1bZWXbV9KNxoo/aP11+/d9/7ADIESpIkSepSVRXQEBfFqa+Pn/c25/XYuHHw1a/GR5MDpgVI/p7oviyuDwEee6znnxsCLF4ch9D//AfefLP1a/Pj556DefPaX1caR6yfPwhMuoP12BXIMhimUk0/bHjlauCwVE7WuxoCJUmSpCGkx60fujM2BVOq4MDZ0DCR9Iwa6urKulgR2b4lQ3crJ6urC3gLXQjx3sTx42H7jC334r2Kb73VGhDnzIGLL+aY54BvfpO3ATbeCHbdNT522SWe6SzNEMWqquKUDZQ3kLO0XfDN4kMICeAnwGrAs1EU/aa7a2wWL0mSpFxpbgJeLCqnV1L3Xj0k0pBOQEM5FTNrO88mbbqj96RReo/OzbZxfXc37vfE3GrYubD8J3ASv+ByTs7JZ3QlQKfN4nMaAkMINwJfBd6PomjLNq/vA1wOlADXR1F0URf3+AawH/AB8Mcoih7p7nMNgZIkScqVQg6ByZlJpj2a4/526RJYuBGs+e8+3aZ6cjXJmbRMA+YlBObiuqZrQxKi6qj1NlEUzxQ+8URc9OaJJ+Bf/4r3H4YQzw6+9RYsXAhRRGOAkkkVvZ4JDCHkLQTuDiwCbm4OgSGEEuBV4IvAXOAZYApxILywwy2OajoWRFF0bQjhriiKDuzucw2BkiRJypVCDoG5UDm9kvqGetJRPBNYsXZ5S7/EjHoyE9jbWcNehsBksgdLT/s7BGby0Ufw9NNxIHziibj4zaJFANSOh8onX+/1nsCuQmBO+wRGUfQY8QxeWzsCc6IoSkVRtAy4DdgviqKXoyj6aofjfeKguKDp2sZcjleSJElSezVTaigfFxeFoaGcmimDpDt7KhWXMoX4ayqVq5WfvbfaarDXXvHGyIcfjmcBX3wRgC2PI2dNEPPRLH594K02z+fSdbmce4AvhRCuBDot6xNCOCaE8GwI4dl5HSv0SJIkSeqVsrFlrTN/02spG9t9MEnOTOZ2UNloU2SlpaRpoSspga23zvnH5CMEZqoP1Ok8axRFn0RRdHQURSdEUXR1F+ddF0XR9lEUbT9+/Ph+GagkSZKknut0X2Iy2dqPr7lsaKbH3R3Nujqnrq613186HT/v6Wdke27BljfNLB8hcC6wYZvnGwDv5GEckiRJkgZSMhlvkGs+IPPj7o5mXZ1TUQGJpriTSMTPe/oZ2Z5rCOzWM8DmIYRNQwjDgW8B2bR5lCRJkqTs1NRAedNexvLy+LmAHIfAEMIM4J/AxBDC3BDC0VEUrQCOBx4CZgF3RFHU/x0QJUmSJBWvsrLW9gq1tTkrsjIYZWhT33+iKJrSyesPAA/k8rMlSZIkSSvLx3JQSZIkSVKeGAIlSZIkqYgYAiVJkiSpiBgCJUmSpB6onlyd7yEMGtX+qAqSIVCSJEnqgeQeyXwPobC1SX6DrH1e0TAESpIkSeo/Jr+CZwiUJEmSlDepFFRWxo8rK+Pnyi1DoCRJkqS8qaqC+vr4cX19/Fy5ZQiUJEmS1KlkEkJofR5C90e25wUi6uognY6vSaehri67+/fkM1yh2p4hUJIkSVKnkkmIotbnUdT9ke15EYGKCkg0pZJEAioqsrt/Tz7DENieIVCSJElSVnLR8qGmBsrL48fl5fFz5VZpvgcgSZIkaXDIxYxaWRnU1sZLN2trs7/OHoS950ygJEmSpEHHJZ69ZwiUJEmSpCJiCJQkSZKkZkXQuNAQKEmSJEnNiqBxoSFQkiRJUuFpblDYp0aEHc7N5rreNC7sqpFhNp85wBscDYGSJEmSCk9zg8I+NSLscG421/WmcWFXjQyz+UxDoCRJkiRlMBBhqQgaFxoCJUmSJA2s3hZfmTYtd2Nq1ty4EOKvZWW5/8wBZgiUJEmSNLCKoPhKITMESpIkSVpZx8Is0LtiLNkWX8nm/s3nZHOu3eQ7ZQiUJEmStLKOhVmgd8VYsi2+ks39m8/J5lxDYKcMgZIkSZIGVhEUXylkpfkegCRJkqQi01x8JYTWIiwaMM4ESpIkSVIRMQRKkiRJUhExBEqSJElSf/tg0161QhwIhkBJkiRJg0JqLFROj5NV5fRKUgsKKFl1NKOmYFshGgIlSZIkDQpVU6C+IU5W9Q31VM3InKzatTgk6lWrw7bX9fgggoaJK7VC7NH1PRlrD1skhqi5z8YQsv3220fPPvtsvochSZIkDRlhWiCqzpwdkjOTTHt02gCPKDeqJ1eT3CMZp6neZqUQCONfITG/knQ6boVYXt6DQqghEJJ0+vPO7hbhuSiKts/0ni0iJEmSJHWrenJ1p+8l90jGwalJV4GxnTZBK5vMVXlcoH7tBOkoTSIkKB9XTu3UbpJVDz+j30ypovwvKerqCq8VostBJUmSJHWrbcjLl5oZUD4ubjJfPq6cmikFlKyapBakqJwKrPlvOK4SxqaorY1bIxYKQ6AkSZKkQaFsAS0zf7VTaykbW0DJqknVjCrqx8WP6xvqYUoBVYRp4p5ASZIkSf0qV8tBm08qlv2JPdWynxH3BEqSJEkqEh33J/Z6T2AvNxBWTq+k/r060glIhATp98qJrs62IkzrZ/e1MExXXA4qSZIkSf2kZkoN5Q3x4/Jx5TCj8PYtGgIlSZIkqZ+UjS2jdnr8uHZqLSwovH2LhkBJkiRJKiKGQEmSJEkqIoZASZIkSSoihkBJkiRJKiKGQEmSJEkqIoZASZIkSSoihkBJkiRJKiKGQEmSJEkqIoZASZIkSSoihkBJkiRJKiKGQEmSJEkqIoZASZIkSeoPqRRUVgLwytVNzwuQIVCSJEnSkNYmm1FZmcNsVlUF9fUAlDc0PS9AhkBJkiRJg18yCSGsfABVE2qpr2sEoL6ukaoJtZnP7XBdt+d0POrqIJ0GoCQift6b+zTL9vxkskc/KkOgJEmSpAHVkteIepa5ms7PeO60ZPx+hqOOStKUAJCmhDoqOz23L0ctFTQ2RazGAFRUxIONop4dzbI93xAoSZIkqZAlk035hdAu93Sbd5rOz+rcNtdUVECiKfkkEnE2y/azenJUvl5DSUU5APXjgJqaAf/ZZsMQKEmSJGlIq6mB8jibUV6ew2xWVga1tQBseVzT8wJUmu8BSJIkSVIuNWezEFoyWlFzJlCSJElSXiVnJvM9hKJiCJQkSZKUV9MenZbvIRQVQ6AkSZIkFRFDoCRJkiQVEUOgJEmSpKGrujrfIyg4hkBJkiRJQ1cPG6kXA0OgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkqR+kVqQonJ6JQCV0ytJLUjleUTKxBAoSZIkqV9UzaiivqEegPqGeqpmVOV5RMqkNN8DkCRJkjSwkjOTTHt0Wk4/Ix2lqZtXR5gWuhgI0Px+MhC6GVL1HvEl6htDoCRJklRkknskSe6R7Pf7Vk6vpL6hnnSUJhESlI8rp3ZqbecXhABRFAfFZEQUdfMByS4CpbLmclBJkiRJ/aJmSg3l48oBKB9XTs2UmjyPSJkYAiVJkiT1i7KxZS0zf7VTaykbW5bnESkTQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVEUOgJEmSJBURQ6AkSZIkFRFDoCRJkiQVkdJ8D6A7IYSNgKuABuDVKIouyvOQJEmSJGnQyulMYAjhxhDC+yGEVzq8vk8IYXYIYU4I4cxubrMF8Mcoio4CKnI2WEmSJEkqArleDvprYJ+2L4QQSoCrgS8Th7opIYSKEMJWIYT7OxyfAV4AvhVC+CvwtxyPV5IkSZKGtJwuB42i6LEQwiYdXt4RmBNFUQoghHAbsF8URRcCX+14jxDCaUB1073uAm7K5ZglSZIkaSgLURTl9gPiEHh/FEVbNj0/ENgniqLvNj0/DNgpiqLjO7l+SyBJvCdwURRFp3Vy3jHAMU1PK4Hafvw2+tvqwIf5HgQDP45cf14u7t+f9xxH/O+xikuh/Pde6Ibaz2kwfD+FMEZ/Dw38ff1dVJwK4b/3wWCo/Zw2j6Jo9Uxv5KMwTMjwWqdJNIqiV4ADu7tpFEXXAdcBhBCui6LomG4uyZtCGd9AjyPXn5eL+/fnPUMIz0ZRtH1/3EuDR6H8917ohtrPaTB8P4UwRn8PDfx9/V1UnArhv/fBYKj9nEII13X2Xj5aRMwFNmzzfAPgnX7+jJp+vl9/K5TxDfQ4cv15ubh/ofyz0uDlv0PZGWo/p8Hw/RTCGP09lN/7qnj471B2htrPqdPvJx/LQUuBV4E9gbeBZ4BDoigq5OWbUp/5f18lSfnm7yJJkPsWETOAfwITQwhzQwhHR1G0AjgeeAiYBdxhAFSR6HRKXpKkAeLvIkm5nwmUJEmSJBWOfOwJlCRJkiTliSFQkiRJkoqIIVCSJEmSioghUCoAIYSvhxB+FUL4Qwhh73yPR5JUfEIIZSGEG0IId+V7LJJyyxAo9VEI4cYQwvshhFc6vL5PCGF2CGFOCOHMru4RRdHvoyj6HnAk8M0cDleSNAT10++iVBRFR+d2pJIKgdVBpT4KIewOLAJubtMPs4S4H+YXgbnE/TCnACXAhR1ucVQURe83XXcZcGsURc8P0PAlSUNAP/8uuiuKogMHauySBl5pvgcgDXZRFD0WQtikw8s7AnOiKEoBhBBuA/aLouhC4Ksd7xFCCMBFwIMGQElST/XH7yJJxcPloFJurA+81eb53KbXOnMCsBdwYAjh+7kcmCSpaPTod1EIYa0Qwi+BbUMIZ+V6cJLyx5lAKTdChtc6XXsdRdEVwBW5G44kqQj19HfRfMD/ESkVAWcCpdyYC2zY5vkGwDt5GoskqTj5u0hSRoZAKTeeATYPIWwaQhgOfAu4L89jkiQVF38XScrIECj1UQhhBvBPYGIIYW4I4egoilYAxwMPAbOAO6Ioqs3nOCVJQ5e/iyT1hC0iJEmSJKmIOBMoSZIkSUXEEChJkiRJRcQQKEmSJElFxBAoSZIkSUXEEChJkiRJRcQQKEmSJElFxBAoSZIkSUXEEChJkiRJRcQQKElSDoQQrgwhPB9C2CGLc8tCCDeEEO4aiLFJkoqbIVCSpH4WQlgV+AxwLPDV7s6PoigVRdHROR+YJElAab4HIElSIQohbABcDVQAJcADwA+jKPq0k/OvBW6OouiJKIoWhxDWBWYCG7U5Zyvgwg6XHhVF0fs5+BYkScrImUBJkjoIIQTgHuD3URRtDmwOjAIu6eKynYAnm65fC1gF+BhobD4hiqKXoyj6aofDAChJGlCGQEmSVvYFYGkURTcBRFHUCJwCHB5CGN3x5BDCJODVpvMAzgF+BtQSzyR2KYSwVgjhl8C2IYSz+ul7kCQpI5eDSpK0skrgubYvRFH0UQjhDWAz4MUO538Z+BNACGETYBfgVGC3pnv9o6sPi6JoPvD9vg9bkqTuORMoSdLKAhB18nomX6IpBAI/Bc6PoigCZhGHQEmSCoYzgZIkrawWOKDtCyGE1YC1gdkdXl8FWCOKondCCNsA+wO7hRCuBkYCLw/IiCVJypIzgZIkrewRYJUQwuEAIYQS4DLgqiiKlnQ49/PA35oeXwxURVG0SRRFmwBb40ygJKnAGAIlSeqgaSnnN4ADQwivAfOBdBRF/y/D6V8G/hRC+AKwahRFj7S5z3vAqiGENQdi3JIkZSPEv+ckSVJnQgi7ADOA/aMoeq7De88DO0VRtDwvg5MkqYcMgZIkSZJURFwOKkmSJElFxBAoSZIkSUXEEChJkiRJRcQQKEmSJElFxBAoSZIkSUXEEChJkiRJRcQQKEmSJElFxBAoSZIkSUXk/wMrYHij/imHzwAAAABJRU5ErkJggg==
"
>
</div>

</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[523]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">siblockfitdata</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2oblock</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2oblock</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2oblock</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2oblock</span><span class="o">.</span><span class="n">x_err</span><span class="p">]</span>
<span class="n">simucinsd2oblockfitdata</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2omucins</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2omucins</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2omucins</span><span class="o">.</span><span class="n">x_err</span><span class="p">]</span>
<span class="n">simucinsh2oblockfitdata</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_h2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_h2omucins</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_h2omucins</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_h2omucins</span><span class="o">.</span><span class="n">x_err</span><span class="p">]</span>
<span class="n">siblockfit</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2oblock</span><span class="o">.</span><span class="n">x</span><span class="p">,</span><span class="n">objective_d2oblock</span><span class="o">.</span><span class="n">generative</span><span class="p">()]</span>
<span class="n">simucinsd2oblockfit</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span><span class="n">objective_d2omucins</span><span class="o">.</span><span class="n">generative</span><span class="p">()]</span>
<span class="n">simucinsh2oblockfit</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_h2omucins</span><span class="o">.</span><span class="n">x</span><span class="p">,</span><span class="n">objective_h2omucins</span><span class="o">.</span><span class="n">generative</span><span class="p">()]</span>

<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;siblockfitdata.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">siblockfitdata</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;simucinsd2oblockfitdata.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">simucinsd2oblockfitdata</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;simucinsh2oblockfitdata.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">simucinsh2oblockfitdata</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;siblockfit.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">siblockfit</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;simucinsd2oblockfit.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">simucinsd2oblockfit</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;simucinsh2oblockfit.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">simucinsh2oblockfit</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>

<span class="kn">import</span> <span class="nn">pickle</span>
<span class="c1"># save</span>
<span class="k">with</span> <span class="nb">open</span><span class="p">(</span><span class="s1">&#39;global_objective.pkl&#39;</span><span class="p">,</span> <span class="s1">&#39;wb+&#39;</span><span class="p">)</span> <span class="k">as</span> <span class="n">f</span><span class="p">:</span>
    <span class="n">pickle</span><span class="o">.</span><span class="n">dump</span><span class="p">(</span><span class="n">global_objective</span><span class="p">,</span> <span class="n">f</span><span class="p">)</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[723]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#hPS</span>

<span class="n">hPS_slab1</span> <span class="o">=</span> <span class="n">hps</span><span class="p">(</span><span class="mf">32.9551</span><span class="p">,</span> <span class="mf">6.44189</span><span class="p">)</span> 
<span class="n">hPS_slab1</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">20</span><span class="p">,</span> <span class="mi">50</span><span class="p">))</span>
<span class="n">hPS_slab1</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;hPS thickness&#39;</span>
<span class="n">hPS_slab1</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">4</span><span class="p">,</span> <span class="mi">8</span><span class="p">))</span>
<span class="n">hPS_slab1</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;hPS roughness&#39;</span>
<span class="n">hPS_slab1</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.0</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.0001</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">))</span>
<span class="n">hPS_slab1</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;hPS solvation&#39;</span>

<span class="n">solv_roughness3</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">22.7224</span><span class="p">,</span> <span class="s1">&#39;hPS/Melinex roughness&#39;</span><span class="p">)</span>
<span class="n">solv_roughness3</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">29</span><span class="p">))</span>

<span class="n">s_hPS</span> <span class="o">=</span> <span class="n">air</span> <span class="o">|</span> <span class="n">hPS_slab1</span> <span class="o">|</span> <span class="n">Melinex</span> <span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness3</span><span class="p">)</span>

<span class="n">model_hPS</span> <span class="o">=</span> <span class="n">ReflectModel</span><span class="p">(</span><span class="n">s_hPS</span><span class="p">,</span><span class="n">scale</span><span class="o">=</span><span class="mf">1.01611</span><span class="p">,</span> <span class="n">bkg</span><span class="o">=</span><span class="mf">1.18178e-06</span><span class="p">,</span> <span class="n">dq_type</span> <span class="o">=</span><span class="s1">&#39;pointwise&#39;</span><span class="p">)</span>
<span class="n">model_hPS</span><span class="o">.</span><span class="n">scale</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.9</span><span class="p">,</span> <span class="mf">1.4</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="n">model_hPS</span><span class="o">.</span><span class="n">bkg</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">4e-9</span><span class="p">,</span> <span class="mf">6e-6</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="c1">#model_hPS.dq.setp(bounds=(4, 6), vary=True)</span>

<span class="n">objective_hPS</span> <span class="o">=</span> <span class="n">Objective</span><span class="p">(</span><span class="n">model_hPS</span><span class="p">,</span> <span class="n">data_hPS</span><span class="p">,</span><span class="n">use_weights</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>

<span class="c1">#Confined mucins 1 Bar</span>

<span class="n">sio2_slab</span> <span class="o">=</span> <span class="n">sio2</span><span class="p">(</span><span class="mf">13.8936</span><span class="p">,</span> <span class="mf">2.19514</span><span class="p">)</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">13.8</span><span class="p">,</span> <span class="mi">14</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;sio2 thickness&#39;</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">1.9</span><span class="p">,</span> <span class="mi">8</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;sio2 roughness&#39;</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.0</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;sio2 solvation&#39;</span>

<span class="n">silanes_slab</span> <span class="o">=</span> <span class="n">silanes</span><span class="p">(</span><span class="mf">23.4608</span> <span class="p">,</span><span class="mf">10.8417</span><span class="p">)</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">22</span><span class="p">,</span> <span class="mi">24</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;silanes thickness&#39;</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">12</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;silanes/sio2 roughness&#39;</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.188672</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.01</span><span class="p">,</span> <span class="mf">0.3</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;silanes solvation&#39;</span>

<span class="n">mucinsnopocket_slab</span> <span class="o">=</span> <span class="n">mucinsnopocket</span><span class="p">(</span><span class="mf">23.9481</span><span class="p">,</span> <span class="mf">11.7125</span><span class="p">)</span>
<span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">15</span><span class="p">,</span> <span class="mi">155</span><span class="p">))</span>
<span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins no pocket thickness 1bar&#39;</span>
<span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>
<span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;mucins/silanes no pocket r 1bar&#39;</span>
<span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.306441</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.01</span><span class="p">,</span> <span class="mf">1.0</span><span class="p">))</span>
<span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins no pocket solvation&#39;</span>

<span class="n">mucinspocket_slab</span> <span class="o">=</span> <span class="n">mucinspocket</span><span class="p">(</span><span class="mf">93.4844</span> <span class="p">,</span> <span class="mf">10.9347</span><span class="p">)</span>
<span class="n">mucinspocket_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">140</span><span class="p">))</span>
<span class="n">mucinspocket_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins pocket thickness&#39;</span>
<span class="n">mucinspocket_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>
<span class="n">mucinspocket_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">constraint</span> <span class="o">=</span> <span class="n">mucinsnopocket_slab</span><span class="o">.</span><span class="n">rough</span>
<span class="c1">#mucinspocket_slab.rough.name = name=&#39;mucins/silanes pocket roughness&#39;</span>
<span class="n">mucinspocket_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.325543</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.01</span><span class="p">,</span> <span class="mf">1.0</span><span class="p">))</span>
<span class="n">mucinspocket_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins pocket solvation&#39;</span>


<span class="n">hPS_slab</span> <span class="o">=</span> <span class="n">hps</span><span class="p">(</span><span class="mf">32.9551</span><span class="p">,</span> <span class="mf">21.9099</span><span class="p">)</span> 
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">20</span><span class="p">,</span> <span class="mi">80</span><span class="p">))</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;hPS thickness&#39;</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">constraint</span> <span class="o">=</span> <span class="n">hPS_slab1</span><span class="o">.</span><span class="n">thick</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">2</span><span class="p">,</span> <span class="mi">60</span><span class="p">))</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;hPS roughness&#39;</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.0</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">))</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;hPS solvation&#39;</span>

<span class="n">solv_roughness4</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">12.3995</span><span class="p">,</span> <span class="s1">&#39;d2o/hPS/mucins roughness&#39;</span><span class="p">)</span>
<span class="n">solv_roughness4</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">3</span><span class="p">,</span> <span class="mi">20</span><span class="p">))</span>

<span class="n">solv_roughness5</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">22.6623</span><span class="p">,</span> <span class="s1">&#39;Melinex/d2o/hPS roughness&#39;</span><span class="p">)</span>
<span class="n">solv_roughness5</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">17</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>

<span class="n">scale1</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">0.6039300181837898</span><span class="p">,</span> <span class="s1">&#39;Scale No pockets&#39;</span><span class="p">)</span>
<span class="n">scale1</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.6</span><span class="p">,</span><span class="mf">0.72</span><span class="p">))</span>
<span class="n">scale2</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">0.007841735744284222</span><span class="p">,</span> <span class="s1">&#39;Scale pockets&#39;</span><span class="p">)</span>
<span class="n">scale2</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">1</span><span class="p">))</span>

<span class="n">s_d2omucinsconfined1</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span> <span class="n">mucinsnopocket_slab</span> <span class="o">|</span> <span class="n">hPS_slab</span> <span class="o">|</span> <span class="n">nopockets</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness4</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined2</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span> <span class="n">mucinspocket_slab</span> <span class="o">|</span> <span class="n">pockets</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness5</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined1</span><span class="o">.</span><span class="n">solvent</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">6.36</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined2</span><span class="o">.</span><span class="n">solvent</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">6.36</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined</span> <span class="o">=</span> <span class="p">[</span><span class="n">s_d2omucinsconfined1</span><span class="p">,</span><span class="n">s_d2omucinsconfined2</span><span class="p">]</span>
<span class="n">scale</span> <span class="o">=</span> <span class="p">[</span><span class="n">scale1</span><span class="p">,</span><span class="n">scale2</span><span class="p">]</span>

<span class="n">model_d2omucinsconfined</span> <span class="o">=</span> <span class="n">MixedReflectModel</span><span class="p">(</span><span class="n">s_d2omucinsconfined</span><span class="p">,</span> <span class="n">scale</span><span class="p">,</span> <span class="n">dq_type</span> <span class="o">=</span><span class="s1">&#39;pointwise&#39;</span><span class="p">)</span>
<span class="n">model_d2omucinsconfined</span><span class="o">.</span><span class="n">bkg</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">4.5921055645868026e-07</span><span class="p">,</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">1e-7</span><span class="p">,</span> <span class="mf">9e-5</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">)</span>
<span class="c1">#model_d2omucinsconfined.dq.setp(8.027801449301215,bounds=(1, 12), vary=False)</span>

<span class="n">objective_d2omucinsconfined</span> <span class="o">=</span> <span class="n">Objective</span><span class="p">(</span><span class="n">model_d2omucinsconfined</span><span class="p">,</span> <span class="n">data_d2omucinsconfined</span><span class="p">,</span><span class="n">use_weights</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="n">global_objective</span> <span class="o">=</span> <span class="n">GlobalObjective</span><span class="p">([</span><span class="n">objective_d2omucinsconfined</span> <span class="p">])</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[726]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">fitter</span> <span class="o">=</span> <span class="n">CurveFitter</span><span class="p">(</span><span class="n">global_objective</span><span class="p">)</span>
<span class="n">fitter</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="s1">&#39;least_squares&#39;</span><span class="p">);</span>
<span class="n">fitter</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="s1">&#39;differential_evolution&#39;</span><span class="p">);</span>

<span class="c1">#fitter = CurveFitter(global_objective, nwalkers=200)</span>
<span class="c1">#np.random.seed(6)</span>
<span class="c1">#fitter.initialise(&#39;jitter&#39;)</span>
<span class="c1">#fitter.reset()</span>

<span class="c1">#fitter.sample(400, random_state=1,pool=2);</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[720]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#fitter.sampler.reset()</span>
<span class="c1">#fitter.sample(10, nthin=200, random_state=1,pool=18);</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="application/vnd.jupyter.stderr">
<pre>100%|██████████████████████████████████████████████████████████████████████████████| 2000/2000 [14:21&lt;00:00,  2.32it/s]
</pre>
</div>
</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[721]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#global_objective.corner();</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>




<div class="jp-RenderedImage jp-OutputArea-output ">
<img src="data:image/png;base64,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
"
>
</div>

</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[725]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#rep = objective_processing.report()</span>
<span class="c1">#rep.process_objective(objective_d2oblock)</span>
<span class="c1">#fig, ax = objective_processing.plot_reports(rÑep, refl_mode=&#39;rq4&#39;)</span>
<span class="c1">#ax[2].set_xscale(&#39;log&#39;)</span>
<span class="c1">#ax[0].get_legend().remove()</span>

<span class="c1">#print(objective_d2oblock)</span>
<span class="c1">#process_chain(global_objective, fitter.chain);</span>
<span class="nb">print</span><span class="p">(</span><span class="n">global_objective</span><span class="p">)</span>
<span class="c1"># the data</span>
<span class="n">plt</span><span class="o">.</span><span class="n">rcParams</span><span class="p">[</span><span class="s1">&#39;figure.figsize&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="p">[</span><span class="mi">15</span><span class="p">,</span> <span class="mi">12</span><span class="p">]</span>

<span class="n">plt</span><span class="o">.</span><span class="n">errorbar</span><span class="p">(</span><span class="n">data_hPS</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_hPS</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_hPS</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_hPS</span><span class="o">.</span><span class="n">x_err</span><span class="p">,</span>
<span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular</span><span class="si">{hPS}</span><span class="s1">$&#39;</span><span class="p">,</span> <span class="n">ms</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">marker</span><span class="o">=</span><span class="s1">&#39;o&#39;</span><span class="p">,</span> <span class="n">lw</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">elinewidth</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;r&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">errorbar</span><span class="p">(</span><span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">x_err</span><span class="p">,</span>
<span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{P=1\ bar}$&#39;</span><span class="p">,</span> <span class="n">ms</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">marker</span><span class="o">=</span><span class="s1">&#39;o&#39;</span><span class="p">,</span><span class="n">lw</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">elinewidth</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;b&#39;</span><span class="p">)</span>

<span class="c1"># the median of the posterior</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">data_hPS</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">objective_hPS</span><span class="o">.</span><span class="n">generative</span><span class="p">(),</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;r&#39;</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">20</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">objective_d2omucinsconfined</span><span class="o">.</span><span class="n">generative</span><span class="p">(),</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;r&#39;</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">20</span><span class="p">,</span> <span class="n">linewidth</span><span class="o">=</span><span class="mi">4</span><span class="p">)</span>

<span class="c1"># plot the spread of the fits for the different datasets</span>
<span class="c1">#gen = objective_d2oblock.pgen(500)</span>
<span class="c1">#save_parsfucins = np.copy(objective_d2omucins.parameters)</span>
<span class="c1">#for i in range(100):</span>
<span class="c1">#    objective_d2omucins.setp(next(gen))</span>
<span class="c1">#    plt.plot(data_d2omucins.x, objective_d2omucins.generative(),</span>
<span class="c1">#    color=&#39;k&#39;, alpha=0.02, zorder=10)</span>
<span class="c1"># put back the saved parameters</span>
<span class="c1">#objective_d2omucins.setp(save_parsfucins)</span>
<span class="n">ax</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">gca</span><span class="p">()</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xscale</span><span class="p">(</span><span class="s2">&quot;log&quot;</span><span class="p">,</span> <span class="n">nonpositive</span><span class="o">=</span><span class="s1">&#39;clip&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_yscale</span><span class="p">(</span><span class="s2">&quot;log&quot;</span><span class="p">,</span> <span class="n">nonpositive</span><span class="o">=</span><span class="s1">&#39;clip&#39;</span><span class="p">)</span>
<span class="c1"># ax.text(-0.04, 1e-11, &#39;a)&#39;)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">legend</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">yscale</span><span class="p">(</span><span class="s1">&#39;log&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylabel</span><span class="p">(</span><span class="s1">&#39;Reflectivity&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlabel</span><span class="p">(</span><span class="s1">&#39;Q /$\AA^{-1}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylim</span><span class="p">(</span><span class="mf">0.5</span><span class="o">*</span><span class="mf">1e-6</span><span class="p">,</span> <span class="mi">2</span><span class="p">);</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlim</span><span class="p">(</span><span class="mf">0.004</span><span class="p">,</span> <span class="mf">0.2</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;Confined1Bar.pdf&#39;</span><span class="p">)</span>


<span class="c1">#np.savetxt(&#39;fitmucinsd2o.txt&#39;,save_parsfucins)</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="text/plain">
<pre>_______________________________________________________________________________

--Global Objective--
________________________________________________________________________________
Objective - 2289394067144
Dataset = d2oSilanesmucinsconfined
datapoints = 181
chi2 = 339.31234334861193
Weighted = True
Transform = None
________________________________________________________________________________
Parameters:       &#39;&#39;       
[Parameters(data=[Parameters(data=[Parameter(value=0.6039300181837898, name=&#39;Scale No pockets&#39;, vary=False, bounds=Interval(lb=0.6, ub=0.72), constraint=None), Parameter(value=0.007841735744284222, name=&#39;Scale pockets&#39;, vary=False, bounds=Interval(lb=0.001, ub=1.0), constraint=None)], name=&#39;scale factors&#39;), Parameter(value=4.5921055645868025e-07, name=&#39;bkg&#39;, vary=False, bounds=Interval(lb=1e-07, ub=9e-05), constraint=None), Parameter(value=5.0, name=&#39;dq - resolution&#39;, vary=False, bounds=Interval(lb=-np.inf, ub=np.inf), constraint=None)], name=&#39;instrument parameters&#39;)]
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;mucins no pocket thickness 1bar&#39;, value=23.9481 (fixed)  , bounds=[15.0, 155.0]&gt;
&lt;Parameter: &#39;sld mucins&#39;  , value=5.55583 (fixed)  , bounds=[5.2, 6.2]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;mucins/silanes no pocket r 1bar&#39;, value=11.7125 (fixed)  , bounds=[1.0, 30.0]&gt;
&lt;Parameter:&#39;mucins no pocket solvation&#39;, value=0.306441 (fixed)  , bounds=[0.01, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;hPS thickness&#39;, value=32.9551          , bounds=[20.0, 80.0], constraint=&lt;Parameter:&#39;hPS thickness&#39;, value=32.9551 (fixed)  , bounds=[20.0, 50.0]&gt;&gt;
&lt;Parameter:   &#39;sld hps&#39;   , value=1.412 (fixed)  , bounds=[1.2, 2.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;hPS roughness&#39;, value=21.9099 (fixed)  , bounds=[2.0, 60.0]&gt;
&lt;Parameter:&#39;hPS solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;no pocket SLD 1&#39;, value=2.61373 (fixed)  , bounds=[2.2, 3.0]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;d2o/hPS/mucins roughness&#39;, value=12.3995 (fixed)  , bounds=[3.0, 20.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:   &#39; - sld&#39;    , value=6.36 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;mucins pocket thickness&#39;, value=93.4844 (fixed)  , bounds=[1.0, 140.0]&gt;
&lt;Parameter: &#39;sld mucins&#39;  , value=5.55583 (fixed)  , bounds=[5.2, 6.2]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=11.7125          , bounds=[0.001, 30.0], constraint=&lt;Parameter:&#39;mucins/silanes no pocket r 1bar&#39;, value=11.7125 (fixed)  , bounds=[1.0, 30.0]&gt;&gt;
&lt;Parameter:&#39;mucins pocket solvation&#39;, value=0.325543 (fixed)  , bounds=[0.01, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;pocket SLD 1&#39; , value=4.94199 (fixed)  , bounds=[4.0, 5.0]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;Melinex/d2o/hPS roughness&#39;, value=22.6623 (fixed)  , bounds=[17.0, 30.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:   &#39; - sld&#39;    , value=6.36 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;


</pre>
</div>
</div>

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>




<div class="jp-RenderedImage jp-OutputArea-output ">
<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA4EAAALCCAYAAABgEda6AAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjMuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8vihELAAAACXBIWXMAAAsTAAALEwEAmpwYAACovElEQVR4nOzdd3hU1dbH8e9OQugiVaSIhJYCiBIFVBQrisTe4Nq9oteCXq++djPYu4Ide1fsBDtW7AZspNACCqJSlEiHZPb7x2aYmcwkhJDkTGZ+n+c5d+bUWSHIzcreey1jrUVEREREREQSQ5LXAYiIiIiIiEj9URIoIiIiIiKSQJQEioiIiIiIJBAlgSIiIiIiIglESaCIiIiIiEgCURIoIiIiIiKSQFK8DqA2GWNygJyWLVue1bt3b6/DERERERER8cT06dOXWWvbRztn4rFPYHZ2ts3Pz/c6DBEREREREU8YY6Zba7OjndN0UBERERERkQSiJFBERERERCSBKAkUERERERFJIHFVGEZERERERBqGjRs3smjRItatW+d1KA1akyZN6NKlC40aNar2PUoCRURERESk3i1atIiWLVuy8847Y4zxOpwGyVrL8uXLWbRoEd27d6/2fZoOKiIiIiIi9W7dunW0bdtWCeA2MMbQtm3brR5NVRIoIiIiIiKeUAK47WryZ6gkUEREREREJIEoCRQRERERkYbD5/M6ggZPSaCIiIiIiDQc48bV2qMWLFhA3759o55LTk5mwIAB9O3bl+OOO441a9YAcOONN5KVlUX//v0ZMGAA33zzTa3FU19UHVRERERERKSCpk2b8sMPPwDwr3/9i4ceeoghQ4YwZcoUZsyYQePGjVm2bBkbNmzwNtAa0EigiIiIiIjEvpISyMpy77Oy3H4tKC8v56yzziIrK4uDDz6YtWvXRlwzdOhQ5s6dy++//067du1o3LgxAO3ataNTp061Ekd9UhIoIiIiIiKxLycHiovd++Jit18L5syZw3nnnUdBQQHbb789r776atj5srIy3nnnHfr168fBBx/MwoUL6d27N+eeey6ffvpprcRQ35QEioiIiIhIbPL5wBi3FRaC3++O+/1uP3BuG4rFdO/enQEDBgAwcOBAFixYAMDatWsZMGAA2dnZ7LTTTpx55pm0aNGC6dOnM3HiRNq3b88JJ5zAk08+uS1foSe0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
"
>
</div>

</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[727]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">sihpsfitdata</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_hPS</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_hPS</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_hPS</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_hPS</span><span class="o">.</span><span class="n">x_err</span><span class="p">]</span>
<span class="n">siconfimucinsh2oblockfitdata</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">x_err</span><span class="p">]</span>
<span class="n">sihpsfit</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_hPS</span><span class="o">.</span><span class="n">x</span><span class="p">,</span><span class="n">objective_hPS</span><span class="o">.</span><span class="n">generative</span><span class="p">()]</span>
<span class="n">siconfimucinsh2oblockfit</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2omucinsconfined</span><span class="o">.</span><span class="n">x</span><span class="p">,</span><span class="n">objective_d2omucinsconfined</span><span class="o">.</span><span class="n">generative</span><span class="p">()]</span>


<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sihpsfitdata.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sihpsfitdata</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;siconfi_P1.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">siconfimucinsh2oblockfitdata</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>

<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sihpsfit.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sihpsfit</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;siconfifit_P1.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">siconfimucinsh2oblockfit</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[810]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#Confined mucins 2 Bar</span>

<span class="n">sio2_slab</span> <span class="o">=</span> <span class="n">sio2</span><span class="p">(</span><span class="mf">13.8936</span><span class="p">,</span> <span class="mf">2.19514</span><span class="p">)</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">13.8</span><span class="p">,</span> <span class="mi">14</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;sio2 thickness&#39;</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">1.9</span><span class="p">,</span> <span class="mi">8</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;sio2 roughness&#39;</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.0</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">))</span>
<span class="n">sio2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;sio2 solvation&#39;</span>

<span class="n">silanes_slab</span> <span class="o">=</span> <span class="n">silanes</span><span class="p">(</span><span class="mf">23.4608</span> <span class="p">,</span><span class="mf">10.8417</span><span class="p">)</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">22</span><span class="p">,</span> <span class="mi">24</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;silanes thickness&#39;</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">12</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;silanes/sio2 roughness&#39;</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.188672</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.01</span><span class="p">,</span> <span class="mf">0.3</span><span class="p">))</span>
<span class="n">silanes_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;silanes solvation&#39;</span>

<span class="n">mucinsnopocket2_slab</span> <span class="o">=</span> <span class="n">mucinsnopocket2</span><span class="p">(</span><span class="mf">23.2952</span><span class="p">,</span><span class="mf">11.7973</span><span class="p">)</span>
<span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span><span class="mf">23.9481</span><span class="p">))</span>
<span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins no pocket thickness 2 bar&#39;</span>
<span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>
<span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;mucins/silanes no pocket r 2bar&#39;</span>
<span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.105266</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">0.306441</span><span class="p">))</span>
<span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins no pocket solvation 2 bar&#39;</span>

<span class="n">mucinspocket2_slab</span> <span class="o">=</span> <span class="n">mucinspocket2</span><span class="p">(</span><span class="mf">48.6421</span><span class="p">,</span> <span class="mf">11.7973</span><span class="p">)</span>
<span class="n">mucinspocket2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">8</span><span class="p">,</span> <span class="mf">93.4844</span><span class="p">))</span>
<span class="n">mucinspocket2_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins pocket thickness 2 bar&#39;</span>
<span class="n">mucinspocket2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>
<span class="n">mucinspocket2_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">constraint</span> <span class="o">=</span> <span class="n">mucinsnopocket2_slab</span><span class="o">.</span><span class="n">rough</span>
<span class="c1">#mucinspocket2_slab.rough.name = name=&#39;mucins/silanes pocket r 2bar&#39;</span>
<span class="n">mucinspocket2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.154216</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">0.325543</span><span class="p">))</span>
<span class="n">mucinspocket2_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;mucins pocket solvation 2 bar&#39;</span>


<span class="n">hPS_slab</span> <span class="o">=</span> <span class="n">hps</span><span class="p">(</span><span class="mf">32.9551</span><span class="p">,</span> <span class="mf">21.9099</span><span class="p">)</span> 
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">20</span><span class="p">,</span> <span class="mi">80</span><span class="p">))</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">thick</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;hPS thickness&#39;</span>
<span class="c1">#hPS_slab.thick.constraint = hPS_slab1.thick</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mi">2</span><span class="p">,</span> <span class="mi">60</span><span class="p">))</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">rough</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="n">name</span><span class="o">=</span><span class="s1">&#39;hPS roughness&#39;</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">0.0279232</span><span class="p">,</span> <span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.0</span><span class="p">,</span> <span class="mf">0.4</span><span class="p">))</span>
<span class="n">hPS_slab</span><span class="o">.</span><span class="n">vfsolv</span><span class="o">.</span><span class="n">name</span> <span class="o">=</span> <span class="s1">&#39;hPS solvation&#39;</span>

<span class="n">scale3</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">0.5847681737764489</span><span class="p">,</span> <span class="s1">&#39;Scale No pockets&#39;</span><span class="p">)</span>
<span class="n">scale3</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.1</span><span class="p">,</span> <span class="mi">1</span><span class="p">))</span>
<span class="n">scale4</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">0.002056453294024637</span><span class="p">,</span> <span class="s1">&#39;Scale pockets&#39;</span><span class="p">)</span>
<span class="n">scale4</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.00001</span><span class="p">,</span> <span class="mi">1</span><span class="p">))</span>

<span class="n">solv_roughness6</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">13.5911</span><span class="p">,</span> <span class="s1">&#39;d2o/hPS/mucins r 2bar&#39;</span><span class="p">)</span>
<span class="n">solv_roughness6</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mf">22.0</span><span class="p">))</span>

<span class="n">solv_roughness8</span> <span class="o">=</span> <span class="n">Parameter</span><span class="p">(</span><span class="mf">15.1441</span><span class="p">,</span> <span class="s1">&#39;Melinex/d2o/hPS r 2bar&#39;</span><span class="p">)</span>
<span class="n">solv_roughness8</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="n">vary</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">0.001</span><span class="p">,</span> <span class="mi">30</span><span class="p">))</span>

<span class="n">s_d2omucinsconfined3</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span> <span class="n">mucinsnopocket2_slab</span> <span class="o">|</span> <span class="n">hPS_slab</span> <span class="o">|</span> <span class="n">nopockets2</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness6</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined4</span> <span class="o">=</span> <span class="n">si</span> <span class="o">|</span> <span class="n">sio2_slab</span> <span class="o">|</span> <span class="n">silanes_slab</span> <span class="o">|</span>  <span class="n">mucinspocket2_slab</span> <span class="o">|</span> <span class="n">pockets2</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">solv_roughness8</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined3</span><span class="o">.</span><span class="n">solvent</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">6.36</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined4</span><span class="o">.</span><span class="n">solvent</span> <span class="o">=</span> <span class="n">SLD</span><span class="p">(</span><span class="mf">6.36</span> <span class="o">+</span> <span class="mi">0</span><span class="n">j</span><span class="p">)</span>
<span class="n">s_d2omucinsconfined2a</span> <span class="o">=</span> <span class="p">[</span><span class="n">s_d2omucinsconfined3</span><span class="p">,</span><span class="n">s_d2omucinsconfined4</span><span class="p">]</span>
<span class="n">scalea</span> <span class="o">=</span> <span class="p">[</span><span class="n">scale3</span><span class="p">,</span><span class="n">scale4</span><span class="p">]</span>

<span class="n">model_d2omucinsconfined2</span> <span class="o">=</span> <span class="n">MixedReflectModel</span><span class="p">(</span><span class="n">s_d2omucinsconfined2a</span><span class="p">,</span> <span class="n">scalea</span><span class="p">,</span> <span class="n">dq_type</span> <span class="o">=</span><span class="s1">&#39;pointwise&#39;</span><span class="p">)</span>

<span class="n">model_d2omucinsconfined2</span><span class="o">.</span><span class="n">bkg</span><span class="o">.</span><span class="n">setp</span><span class="p">(</span><span class="mf">5.979137923565975e-08</span><span class="p">,</span><span class="n">bounds</span><span class="o">=</span><span class="p">(</span><span class="mf">1e-9</span><span class="p">,</span> <span class="mf">2e-5</span><span class="p">),</span> <span class="n">vary</span><span class="o">=</span><span class="kc">True</span><span class="p">)</span>
<span class="c1">#model_d2omucinsconfined2.dq.setp(8.396173635629658,bounds=(1, 12), vary=False)</span>

<span class="n">objective_d2omucinsconfined2</span> <span class="o">=</span> <span class="n">Objective</span><span class="p">(</span><span class="n">model_d2omucinsconfined2</span><span class="p">,</span> <span class="n">data_d2omucinsconfined2</span><span class="p">,</span><span class="n">use_weights</span> <span class="o">=</span> <span class="kc">True</span><span class="p">)</span>
<span class="n">global_objective2</span> <span class="o">=</span> <span class="n">GlobalObjective</span><span class="p">([</span><span class="n">objective_d2omucinsconfined2</span><span class="p">])</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[805]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">fitter2</span> <span class="o">=</span> <span class="n">CurveFitter</span><span class="p">(</span><span class="n">global_objective2</span><span class="p">)</span>
<span class="n">fitter2</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="s1">&#39;least_squares&#39;</span><span class="p">);</span>
<span class="n">fitter2</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="s1">&#39;differential_evolution&#39;</span><span class="p">);</span>

<span class="c1">#fitter2 = CurveFitter(global_objective2, nwalkers=200)</span>
<span class="c1">#np.random.seed(6)</span>
<span class="c1">#fitter2.initialise(&#39;jitter&#39;)</span>
<span class="c1">#fitter2.reset()</span>

<span class="c1">#fitter2.sample(400, random_state=1,pool=2);</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="application/vnd.jupyter.stderr">
<pre>22it [00:04,  4.80it/s]
</pre>
</div>
</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[771]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#fitter2.sampler.reset()</span>
<span class="c1">#fitter2.sample(5, nthin=200, random_state=1,pool=18);</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="application/vnd.jupyter.stderr">
<pre>100%|██████████████████████████████████████████████████████████████████████████████| 1000/1000 [06:06&lt;00:00,  2.73it/s]
</pre>
</div>
</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[800]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#global_objective2.corner();</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[811]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#rep = objective_processing.report()</span>
<span class="c1">#rep.process_objective(objective_d2oblock)</span>
<span class="c1">#fig, ax = objective_processing.plot_reports(rÑep, refl_mode=&#39;rq4&#39;)</span>
<span class="c1">#ax[2].set_xscale(&#39;log&#39;)</span>
<span class="c1">#ax[0].get_legend().remove()</span>

<span class="c1">#print(objective_d2oblock)</span>

<span class="c1"># the data</span>
<span class="n">plt</span><span class="o">.</span><span class="n">rcParams</span><span class="p">[</span><span class="s1">&#39;figure.figsize&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="p">[</span><span class="mi">15</span><span class="p">,</span> <span class="mi">12</span><span class="p">]</span>


<span class="n">plt</span><span class="o">.</span><span class="n">errorbar</span><span class="p">(</span><span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">x_err</span><span class="p">,</span>
<span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{P=2\ bar}$&#39;</span><span class="p">,</span> <span class="n">ms</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">marker</span><span class="o">=</span><span class="s1">&#39;o&#39;</span><span class="p">,</span><span class="n">lw</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">elinewidth</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;b&#39;</span><span class="p">)</span>

<span class="c1"># the median of the posterior</span>

<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">objective_d2omucinsconfined2</span><span class="o">.</span><span class="n">generative</span><span class="p">(),</span> <span class="n">color</span><span class="o">=</span><span class="s1">&#39;r&#39;</span><span class="p">,</span> <span class="n">zorder</span><span class="o">=</span><span class="mi">20</span><span class="p">,</span> <span class="n">linewidth</span><span class="o">=</span><span class="mi">4</span><span class="p">)</span>

<span class="c1"># plot the spread of the fits for the different datasets</span>
<span class="c1">#gen = objective_d2oblock.pgen(500)</span>
<span class="c1">#save_parsfucins = np.copy(objective_d2omucins.parameters)</span>
<span class="c1">#for i in range(100):</span>
<span class="c1">#    objective_d2omucins.setp(next(gen))</span>
<span class="c1">#    plt.plot(data_d2omucins.x, objective_d2omucins.generative(),</span>
<span class="c1">#    color=&#39;k&#39;, alpha=0.02, zorder=10)</span>
<span class="c1"># put back the saved parameters</span>
<span class="c1">#objective_d2omucins.setp(save_parsfucins)</span>
<span class="n">ax</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">gca</span><span class="p">()</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_xscale</span><span class="p">(</span><span class="s2">&quot;log&quot;</span><span class="p">,</span> <span class="n">nonpositive</span><span class="o">=</span><span class="s1">&#39;clip&#39;</span><span class="p">)</span>
<span class="n">ax</span><span class="o">.</span><span class="n">set_yscale</span><span class="p">(</span><span class="s2">&quot;log&quot;</span><span class="p">,</span> <span class="n">nonpositive</span><span class="o">=</span><span class="s1">&#39;clip&#39;</span><span class="p">)</span>
<span class="c1"># ax.text(-0.04, 1e-11, &#39;a)&#39;)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">legend</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">yscale</span><span class="p">(</span><span class="s1">&#39;log&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylabel</span><span class="p">(</span><span class="s1">&#39;Reflectivity&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlabel</span><span class="p">(</span><span class="s1">&#39;Q /$\AA^{-1}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylim</span><span class="p">(</span><span class="mf">0.1</span><span class="o">*</span><span class="mf">1e-6</span><span class="p">,</span> <span class="mi">2</span><span class="p">);</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlim</span><span class="p">(</span><span class="mf">0.004</span><span class="p">,</span> <span class="mf">0.2</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">savefig</span><span class="p">(</span><span class="s1">&#39;Confined2Bar.pdf&#39;</span><span class="p">)</span>

<span class="nb">print</span><span class="p">(</span><span class="n">global_objective2</span><span class="p">)</span>
<span class="c1">#np.savetxt(&#39;fitmucinsd2o.txt&#39;,save_parsfucins)</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>


<div class="jp-RenderedText jp-OutputArea-output" data-mime-type="text/plain">
<pre>_______________________________________________________________________________

--Global Objective--
________________________________________________________________________________
Objective - 2289412519368
Dataset = d2oSilanesmucinsconfined2
datapoints = 206
chi2 = 470.55607845892257
Weighted = True
Transform = None
________________________________________________________________________________
Parameters:       &#39;&#39;       
[Parameters(data=[Parameters(data=[Parameter(value=0.5847681737764489, name=&#39;Scale No pockets&#39;, vary=False, bounds=Interval(lb=0.1, ub=1.0), constraint=None), Parameter(value=0.002056453294024637, name=&#39;Scale pockets&#39;, vary=False, bounds=Interval(lb=1e-05, ub=1.0), constraint=None)], name=&#39;scale factors&#39;), Parameter(value=5.979137923565975e-08, name=&#39;bkg&#39;, vary=True, bounds=Interval(lb=1e-09, ub=2e-05), constraint=None), Parameter(value=5.0, name=&#39;dq - resolution&#39;, vary=False, bounds=Interval(lb=-np.inf, ub=np.inf), constraint=None)], name=&#39;instrument parameters&#39;)]
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;mucins no pocket thickness 2 bar&#39;, value=23.2952 (fixed)  , bounds=[0.001, 23.9481]&gt;
&lt;Parameter: &#39;sld mucins&#39;  , value=5.55583 (fixed)  , bounds=[5.2, 6.2]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;mucins/silanes no pocket r 2bar&#39;, value=11.7973 (fixed)  , bounds=[0.001, 30.0]&gt;
&lt;Parameter:&#39;mucins no pocket solvation 2 bar&#39;, value=0.105266 (fixed)  , bounds=[0.001, 0.306441]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;hPS thickness&#39;, value=32.9551          , bounds=[20.0, 80.0]&gt;
&lt;Parameter:   &#39;sld hps&#39;   , value=1.412 (fixed)  , bounds=[1.2, 2.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;hPS roughness&#39;, value=21.9099          , bounds=[2.0, 60.0]&gt;
&lt;Parameter:&#39;hPS solvation&#39;, value=0.0279232 (fixed)  , bounds=[0.0, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;no pocket SLD 2&#39;, value=2.60686 (fixed)  , bounds=[2.2, 3.0]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;d2o/hPS/mucins r 2bar&#39;, value=13.5911 (fixed)  , bounds=[0.001, 22.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:   &#39; - sld&#39;    , value=6.36 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
________________________________________________________________________________
Parameters: &#39;Structure - &#39; 
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=2.07 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;sio2 thickness&#39;, value=13.8936 (fixed)  , bounds=[13.8, 14.0]&gt;
&lt;Parameter:   &#39; - sld&#39;    , value=3.47 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;sio2 roughness&#39;, value=2.19514 (fixed)  , bounds=[1.9, 8.0]&gt;
&lt;Parameter:&#39;sio2 solvation&#39;, value=0 (fixed)  , bounds=[0.001, 0.4]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;silanes thickness&#39;, value=23.4608 (fixed)  , bounds=[22.0, 24.0]&gt;
&lt;Parameter: &#39;silanes SLD&#39; , value=0.7 (fixed)  , bounds=[0.6, 0.8]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;silanes/sio2 roughness&#39;, value=10.8417 (fixed)  , bounds=[1.0, 12.0]&gt;
&lt;Parameter:&#39;silanes solvation&#39;, value=0.188672 (fixed)  , bounds=[0.01, 0.3]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:&#39;mucins pocket thickness 2 bar&#39;, value=48.6421 (fixed)  , bounds=[8.0, 93.4844]&gt;
&lt;Parameter: &#39;sld mucins&#39;  , value=5.55583 (fixed)  , bounds=[5.2, 6.2]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:  &#39; - rough&#39;   , value=11.7973          , bounds=[0.001, 30.0], constraint=&lt;Parameter:&#39;mucins/silanes no pocket r 2bar&#39;, value=11.7973 (fixed)  , bounds=[0.001, 30.0]&gt;&gt;
&lt;Parameter:&#39;mucins pocket solvation 2 bar&#39;, value=0.154216 (fixed)  , bounds=[0.001, 0.325543]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:  &#39; - thick&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;pocket SLD 2&#39; , value=5.54975 (fixed)  , bounds=[4.6, 6.36]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:&#39;Melinex/d2o/hPS r 2bar&#39;, value=15.1441 (fixed)  , bounds=[0.001, 30.0]&gt;
&lt;Parameter:&#39; - volfrac solvent&#39;, value=0 (fixed)  , bounds=[0.0, 1.0]&gt;
________________________________________________________________________________
Parameters:       &#39;&#39;       
&lt;Parameter:   &#39; - sld&#39;    , value=6.36 (fixed)  , bounds=[-inf, inf]&gt;
&lt;Parameter:   &#39; - isld&#39;   , value=0 (fixed)  , bounds=[-inf, inf]&gt;


</pre>
</div>
</div>

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>




<div class="jp-RenderedImage jp-OutputArea-output ">
<img src="data:image/png;base64,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
"
>
</div>

</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[807]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">siconfimucinsh2oblockfitdata1</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">x</span><span class="p">,</span> <span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">y</span><span class="p">,</span> <span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">y_err</span><span class="p">,</span><span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">x_err</span><span class="p">]</span>
<span class="n">siconfimucinsh2oblockfit1</span> <span class="o">=</span> <span class="p">[</span><span class="n">data_d2omucinsconfined2</span><span class="o">.</span><span class="n">x</span><span class="p">,</span><span class="n">objective_d2omucinsconfined2</span><span class="o">.</span><span class="n">generative</span><span class="p">()]</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;siconfi_P2.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">siconfimucinsh2oblockfitdata1</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;siconfifit_P2.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">siconfimucinsh2oblockfit1</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell   ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[821]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#SLD</span>

<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="o">*</span><span class="n">s_d2omucinsconfined1</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">(),</span><span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{P=1\ bar}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="o">*</span><span class="n">s_d2omucinsconfined2</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">(),</span><span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{P=2\ bar}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="o">*</span><span class="n">s_d2omucins</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">(),</span><span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{Si\ mucins\ D_2O}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="o">*</span><span class="n">s_h2omucins</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">(),</span><span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{Si\ mucins\ H_2O}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="o">*</span><span class="n">s_d2oblock</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">(),</span><span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular{Si\ block}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">plot</span><span class="p">(</span><span class="o">*</span><span class="n">s_hPS</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">(),</span><span class="n">label</span><span class="o">=</span><span class="s1">&#39;$\mathregular</span><span class="si">{hPS}</span><span class="s1">$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylim</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="mi">7</span><span class="p">);</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlim</span><span class="p">(</span><span class="o">-</span><span class="mi">40</span><span class="p">,</span> <span class="mi">800</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">legend</span><span class="p">()</span>
<span class="n">plt</span><span class="o">.</span><span class="n">ylabel</span><span class="p">(</span><span class="s1">&#39;SLD /$10^{-6} \AA^{-2}$&#39;</span><span class="p">)</span>
<span class="n">plt</span><span class="o">.</span><span class="n">xlabel</span><span class="p">(</span><span class="s1">&#39;distance / $\AA$&#39;</span><span class="p">);</span>
</pre></div>

     </div>
</div>
</div>
</div>

<div class="jp-Cell-outputWrapper">


<div class="jp-OutputArea jp-Cell-outputArea">

<div class="jp-OutputArea-child">

    
    <div class="jp-OutputPrompt jp-OutputArea-prompt"></div>




<div class="jp-RenderedImage jp-OutputArea-output ">
<img src="data:image/png;base64,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ADDU+vXr9fjjj2vjxo168803tW7dOhUXF0uSXn755UmrI1rPSklJ0ebNm7V161b94x//0F//+ld97WtfG/Z6a63e//736/LLL9euXbu0c+dOtba26stf/nJU6hmLsXQKT5L0LmvtVdbaayQtiXFNU1drjdRWNznPajooSWpLyZIU3ry+qjG8V+DMrMkNhfPz0uV1GzaxBwAAwBCVlZXKzc1VUlJ4McTc3FwVFhZKktLT04dcX15eriVLluijH/2oli1bpg996ENat26dzjjjDC1cuFCvvfZa33XLli3re9/dd9+tr371q5KkBx54QCtWrNDKlSt1/fXXD3hWeXm5SktLdfPNN2vp0qW66KKL1NERbq60tbXpXe96l1auXKlly5bpd7/73Yjf24wZM7R27Vr96Ec/Oro93CBPP/20kpOT9eEPf1hSuNP4ve99T/fee6/a29vH9DOcqFFXH5XULGmWMeZmSUckTc6u51NNS5We+N6/KjPUqTP+9f+TClbE9nk9obDVF/6PO82bppqmThVJ8k3SdhS9vG6XFszI0NuVdAoBAAAS1bdf+7a2N2yP6j2X5CzRF07+wojXXHTRRfr617+uRYsW6YILLtDVV1+tc845Z8T37N69W7///e+1du1anXTSSfrNb36jF198UX/+859111136dFHHx32vVu3btWdd96pl156Sbm5uWpoaBhyza5du/Tggw/qZz/7mT7wgQ/o4Ycf1nXXXacnn3xShYWF+stf/iJJampqGvVnMG/ePIVCIdXU1GjmzJkR61mzZs2AY5mZmZo9e7Z2796tFStinBs0tk7hVyT9SVKOJJ+k22Ja0RR1+G9/1N7Dn9bmqi/o8B//N/YPbKqQ3Elq94Rzf7o3XTXN4TmFxpjYP3+Q0gJCIQAAAIZKT0/Xhg0btHbtWuXl5enqq6/WfffdN+J7SkpKtHz5crlcLi1dulTnn3++jDFavny5ysvLR3zv008/rSuvvFK5ubmSpJycnIj3LysrkyStWbOm757Lly/XunXr9IUvfEEvvPCCsrKyxvQ99nYJH330Ud1888267LLL9Pe//73vXKS/z4c7Hgtj2afQSnrUGJNrrZ2ksY9Tz4YN7bIKqdvdpad3zdd1jRVSdnHsHth0SMoqUqs/vLhLqidV1c2dkjxyuSY/FJ5QkKk/bjykutYu5aZP3j6JAAAAGJvROnqx5Ha7de655+rcc8/V8uXLdf/99+vGG28c9vreoaaS5HK5+r52uVwKBAKSJI/HM2CBl87OcINkLGGr//3dbnff8NFFixZpw4YN+utf/6ovfelLuuiii3THHXeMeK+9e/fK7XZrxowZuvzyy3X55ZfryJEjuv3223XRRRdp6dKlevjhhwe8p7m5WRUVFZo/f/6I946W8aw+em/MqpjqQiHVtBSofNoWbc1/QY2dK9S28a+xfWbTwZ5Q2CpJ6vb71NkdDJ+b/EyoEwoyJYluIQAAAAbYsWOHdu3a1ff15s2bNWfOnAnfd+bMmaqpqVF9fb26urr0+OOPS5LOP/98PfTQQ6qvr5ekiMNHh3P48GGlpqbquuuu0+23366NGzeOeH1tba0+9rGP6ZOf/OSAIPrNb35Tn/jEJ/rqaW9v1wMPPCBJCgaD+rd/+zfdeOONSk1NHdf3fKzGMqewVxyixNTQXVOuzsBM1aZvVGe+X+aQS2++slOnnRfDhzYdlOafp7aeTmFts+37BcZn+OjRUHjWwrxJfz4AAAASU2trq2677TY1NjbK4/FowYIFWrt27YTv6/V6dccdd+iUU05RSUmJliwJr5e5dOlSffnLX9Y555wjt9utVatWjTpctddbb72lz3/+83K5XPJ6vfrJT34y5JqOjg6VlZXJ7/fL4/Ho+uuv1+c+9zlJ4S7lF7/4RV166aVavXq1pPDf5o888ohuvfVWfeMb31AoFNI73/lO3XXXXRP+GYzVeEJh5OVyMKqKN7dLStaR1Crddvbt2vjGPm1uyNNp7Q1S6tAxzBMW9EstlVJWkdr8bfK5fDrUEDgaCie0O+WxmZbmU35msrYdplMIAACAo9asWTPsdhCtra1Djs2dO1dbtmzp+7p/oBt87lOf+pQ+9alPDbnHDTfcoBtuuCHiswbf4/bbb+/7/OKLL9bFF1884vcTDAaHPffDH/5Q69atU1NTk3bv3q2PfexjkqTi4mI99thjI943lugUToKa/dWS5siX5dIFC8r0RMazSmpcpsCup+VZeWX0H9h8WJKVsmaprfuA0n3pKq9vk7unQxiPTqEknVCYybYUAAAAcKzhQmq8jadn9KWYVdHDGJNsjHnNGPOGMWarMWb4XR6PI8314f1FcvNyZYyRqzBP3mCqXnw9Rptx9mxH0TunMNWTqn11bZqe5pMkxSkTqrQgQ3tqW9XpH/5fTwAAAABMrjGHQmvtFmPMEmPMF4wxPzDGfL/n89Io1tMl6Txr7UpJZZIuMcacGsX7x0V7c1Cd7jYtKZgnSTptzSpJ0vaKgDTMJpYT0hcKi9XW3dbXKczrWfUzXp3C0oJMBUJWu2uGDgMAAAAAEB9jDoXGmC9I+q3Cw0hfk/TPns8fNMZ8MRrF2LDexODt+Tju5zI2d3rVmnREC3PCofCClYvUlFSrpra5UsPe6D+wqSL8mjlLbYE2pXnTVF7XrhnpPZ3COGxJIR1dgXQbK5ACAAAACWM8cwpvkrTUWuvvf9AY811JWyV9KxoFGWPckjZIWiDpx9baV6Nx33jqCKSpNe2I5mWdJEnKSvGqIb1BhY0LFdr1lFzTo7z/SOMBKTVX8qWqtbtVWb5ctXYFlJuepG61xG346JzpaUrxutmWAgAAAEgg45lTGJJUGOF4Qc+5qLDWBq21ZZKKJJ1sjFk2+BpjzC3GmNeNMa/X1tZG69ExEwhmqN3bojmZ/fZbmZGmpGCa3nrr9eg/8Mg+KadEktTmb5MNhTuE09N6ho/GqVPodhktzs9gBVIAAAAggYwnFH5G0lPGmCeMMWt7Pp6U9FTPuaiy1jZKelbSJRHOrbXWnmitPTEvL7H3vLP+bimYLnk6leo9uvnk8hWLJElvHAxIoahl6rAj5dK0uZKkVn+r/P5wKMxN7V1oJn4LyYZXIG2WjcVcSgAAAADjNp6FZp6UtEjS1yT9TdLfJX1V0mJr7RPRKMYYk2eMye75PEXSBZK2R+Pe8dJRXycjt3zegcHvPSeXqc3bpNrO+VL1W9F7YKA7vNDMtHCnsN3frs4urzwuo+wUr6T47FPYq7QgU82dAR1u6oxfEQAAAAD6jGdOoay1IUmvDD5ujPmwtfaXUainQNL9PfMKXZIestY+HoX7xk31oUpJUlrywB91Tlqy6tJqlde+UNr3vFSwMjoPbKqQbEiaNlf+kF+dwU61dng0e3qqekeNxrVTWJAhSdp2uFmzslPiVgcAAACAsGj1jKKyn6C19k1r7Spr7Qpr7TJr7dejcd94Kj9YLkmalp4x5Fwgx63U7lxVbl8fvQce2Rd+zSlRuz+8P2Jjq0vzctP6RqnGs1O4OD9TxojFZgAAAIAEMZ4tKd4c5uMtSTNjWONxraK2SpKUPz1/yLnieeF1e1472CIF/UPOH5OGnlA4ba5a/eHdPepbjUpy0/rm8cWzU5ie5NGcnFRCIQAAAPrceeedWrp0qVasWKGysjK9+mp4A4LTTz990mqI1rPS09MHfH3ffffpk5/85LDXHzx4UJdddpkWLlyo+fPn69Of/rS6u7ujUstYjadnNFPSv0h6T4SP+uiXNjU0tISDWcmskiHnLj51jUIKan/3POnQxug88Ei55EmW0vPV5m+TJAX8PpXkpvft+BjHTCgpPK+QvQoBAAAgSevXr9fjjz+ujRs36s0339S6detUXFwsSXr55ZcnrY7JfFYva63e//736/LLL9euXbu0c+dOtba26stf/vKk1jGeUPi4pHRr7f5BH+UKrxKKCLo7wmM2i4qHhsLSWXk6klKn1q4F4XmF0dC78qjL1RcKbSgpYTqFUngT+/317WrtCsS1DgAAAMRfZWWlcnNzlZQU3j4tNzdXhYXhEXWDu26SVF5eriVLluijH/2oli1bpg996ENat26dzjjjDC1cuFCvvfZa33XLlh3d3e7uu+/WV7/6VUnSAw88oBUrVmjlypW6/vrrBzyrvLxcpaWluvnmm7V06VJddNFF6ujokCS1tbXpXe96l1auXKlly5bpd7/73YS+96efflrJycn68Ic/LElyu9363ve+p3vvvVft7e0Tuvd4jHmhGWvtTSOcuzY65Uw9oW63/K4uzcwbOnzUGKPmjA4V1pcosPcn8pzz+Yk/sP92FN3hLqUNJmteXpoO7Ql/Ha99CnuVFmRKknZUNWvNnJy41gIAAICwqrvuUtfb0V34P6l0ifL/4z9GvOaiiy7S17/+dS1atEgXXHCBrr76ap1zzjkjvmf37t36/e9/r7Vr1+qkk07Sb37zG7344ov685//rLvuukuPPvrosO/dunWr7rzzTr300kvKzc1VQ0PDkGt27dqlBx98UD/72c/0gQ98QA8//LCuu+46PfnkkyosLNRf/vIXSVJTU9OQ93Z0dKisrKzv64aGBr33ve8dtpY1a9YMOJaZmanZs2dr9+7dWrFixYg/h2gZz5zC00y8W0zHoZDfLb+7Q0let6q//R0deeihAedTC1KVFEzVloOHJX/HxB5m7YBQ2BYIdwqT3amakZGUOMNHC8OhkE3sAQAAkJ6erg0bNmjt2rXKy8vT1Vdfrfvuu2/E95SUlGj58uVyuVxaunSpzj//fBljtHz5cpWXl4/43qefflpXXnmlcnNzJUk5OUObFCUlJX3Bbs2aNX33XL58udatW6cvfOELeuGFF5SVlTXkvSkpKdq8eXPfx9e/fnTtzEcffVQ333yzLrvsMv3973+XtTbiKL7hjsfKeLakuEHSj40xOyU9KelJa21VbMqaQoIeBdydCrW1qeFXv1LK0qWa9oEP9J1efsJ81bzRom2BBSqreFWad+6xP6v5sNTdKk1fIElq6w6HwuLsaTLGJMzw0cKsZGWleJlXCAAAkEBG6+jFktvt1rnnnqtzzz1Xy5cv1/33368bb7xx2Ot7h5pKksvl6vva5XIpEAhPUfJ4PAqFju4V3tkZ3id7LIGr//3dbnff8NFFixZpw4YN+utf/6ovfelLuuiii3THHXeM+fu8/PLLdfnll+vIkSO6/fbbdc011+jhhx8ecE1zc7MqKio0f/78Md93osazef3HrLWrFd6wfpqk+4wx640xdxljzu7ZWxCDmKBPQXen2jdtlgIBde7cKRsM9p2/5OQV8ru6VOlfNPF5hXU7wq95iyWpb/XRkp5//bChnlAYxy0ppHAoXT4rS28eHNpuBwAAgLPs2LFDu3bt6vt68+bNmjNnzoTvO3PmTNXU1Ki+vl5dXV16/PHw9ufnn3++HnroIdXXh9fKjDR8dDiHDx9WamqqrrvuOt1+++3auPHYFov85je/qU984hM6//zz1d7ergceeECSFAwG9W//9m+68cYblZqaekz3Phbj2rxekqy12yVtl/Q9Y0yKpHdIukrSdyWdGN3yjn8mmKyQp0PtPcvq2o4Ode8/oKR54YVnMlKSVJ9Wp5TOBdK+X0/sYXU9/2PKDYfC5q5wKFyQ2xMK+4aPxn8U8IqiLK19fq86/UEle/n3BAAAAKdqbW3VbbfdpsbGRnk8Hi1YsEBr166d8H29Xq/uuOMOnXLKKSopKdGSJUskSUuXLtWXv/xlnXPOOXK73Vq1atWow1V7vfXWW/r85z8vl8slr9ern/zkJ+OqyVqrL37xi7r00ku1evVqSdIjjzyiW2+9Vd/4xjcUCoX0zne+U3fddde47jtR4w6F/VlrOyT9tecDEbiCyQr5mtT+2mtyZ2cr2Nioru1v94VCSeqa1q28imJ1HHxTKZ3NUnLmsT2sdoeUnCWlz5AkVbU0yQZ9mpcXvl9vp1Dxz4RaUZSlQMhqe1WLyoqz410OAAAA4mTNmjXDbgfR2to65NjcuXO1ZcuWvq/7B7rB5z71qU/pU5/61JB73HDDDbrhhhsiPmvwPW6//fa+zy+++GJdfPHFI34/g2u+8cYb+4bC/vCHP9S6devU1NSk3bt362Mf+5iKi4v12GOPjXjPWIvzQMKpzx1Mkct0q2PLFmW9//2S16vOQas65RZPk9t69KadK+2fwP4odTul3EV9K8lUtzbKhpJVkpsmKdwpNCYxOoXLi7IlSW8dbIxrHQAAAMBk+dSnPqUNGzbopz/9qT72sY/Fu5w+hMJYCgbkCSYrvaNLCgaVdsbpSpo/X53bB4bCM1aXSpK22yUTm1dYt7Nv6Kgk1Xc09+1RKIU7hYkQCKXwYjPT03zMKwQAAADi7JhCoTEmzxiTF+1ipppAW6M81qeslg7J61XqqlVKXrJEndvfHnDdaaWL1elu16HABEJhR6PUWi3lLeo71NTZIo+SlZ3qk9QzpzBB/hnAGKPlRVl66xChEAAAAIin8exTaIwxXzXG1Cm80MxOY0ytMWbsa7A6TEN9eMeOaUdalbJ8uVypqUouXaJgbZ0CtbV917lcLjWmH1Gws0iqfktqqx//w+p2hl/7dQpbutuU6knr+9paK1eCdAolacWsLO2sblFHd3D0iwEAAADExHj6Rp+RdIakk6y106210ySdIukMY8xnY1Hc8a66ulqSlFXfpNSTT5IkJS0JDxXt3L5jwLVmulVWR77a5JXKj6FbWNtzv9yFksIBsDPQroyk9L5LbMhKrsQJhcuLshWy0rZKuoUAAABAvIwnFP6LpGustft6D1hr90q6ruccBqlrCHf8vIEOpZ18siQpeUm4kzd4CGlxSa7c1qPX3cc4hLR6i+RNk6bNDX/Z3KWQq015qdP6LuldaCZRrCjKkiTmFQIAAABxNJ5Q6LXW1g0+aK2tleSNXklTR1NzsyTJE+hQUs++KO6sLHkKC9Q1aAXSs9cslSRt1ppjC4VVb0kzl0qu8J5/O6ubZTxtKszM7bvE2sRZaEaSZmYma2Zmkt4iFAIAAABxM55Q2H2M5xyrpbVdUjgUujOP7j2YvKR0yAqkC+cUq9vdqZrO2VL9bqnp0NgfZG04FOYv7zu0papWxgQ1b9rMAZeZBFloptfyWdl6k8VmAAAAgLgZT0RYaYxpjvDRImn5qO92oPb2LkmSyxOU8Xj6jicvWaLuffsUam/vO+Zyu9SccUSejvDG8+PqFjbul7qaB4TC7TWVkqTCzOl9xxJpS4peK4qytKe2Va1dgXiXAgAAADjSmEOhtdZtrc2M8JEh6d9jWONxq6szvKqmJ9kz4HjSksWSteras3fAcU+elNMxQ02+3PGFwso3w68FK/oO7a4PL3KTk5zTdyzR5hRK4VBorfRmRWO8SwEAAAAcKVqDCVl9NIJAt5UkJaUnDTjumzNHkuQ/WDHg+LyFefJYr55OOjMcCq0d24MOvS65vNKM8LxEa60qGmskSdOS+i80Y2USaPVRSVpVHK5v44Ejca4EAAAA8XLnnXdq6dKlWrFihcrKyvTqq69Kkk4//fQh15aXl2vZsmUR75Oenh7x+GhGuqcTRCsUJlbSSBABv0sm1Kmk7OwBx72ziiRJ3QcPDji+Zll44/lNnSVS80GpYWAncVgVr0mFZZI3WZJU29KljlB4kZtpyf1CYQIOH81K9WrRzHRt2E8oBAAAcKL169fr8ccf18aNG/Xmm29q3bp1Ki4uliS9/PLLca7OGTyjXzImY2xpOUso4JE72CFvdtaA4+70NLmzs+U/OHAxmQVzZusv7m1qbJsu+STte06aPn/khwS6pcObpBNv6jv0dlWLjDs8X3FAKEzA4aOStGbONP3lzUqFQlauBOtkAgAAOMULD+1UXUVrVO+ZW5yusz6waMRrKisrlZubq6Sk8Oi63Nyjq+enp6ertXVoTYFAQDfccIM2bdqkRYsW6YEHHlBqauqAa7773e/q3nvvlSR99KMf1Wc+8xlJ0gMPPKC7775bxhitWLFCv/rVrwa8b+/evbriiiu0du1anXTSSeP+no9HY+4UGmNaBi8w0/shqTCGNR63bNAjT6BT7ozMIee8RUXyD+oUutwutWU0KqszW/70QmnXutEfUv2WFOiUik/uO7TtcLOMu00+l0+pnn7/4wgl3vBRSVo9e5qaOwPaUxvd/ycEAACAxHfRRRepoqJCixYt0q233qrnnntu1Pfs2LFDt9xyi958801lZmbqnnvuGXB+w4YN+uUvf6lXX31Vr7zyin72s59p06ZN2rp1q+688049/fTTeuONN/T9739/yH2vuOIK/fKXv3RMIJTG0SnsWVAG4xFMktffIXdW5FDY+fa2IceT843Sd+ZpW847tHL3H6TOZil56Pv77Hsh/Dr71L5D2yqblZ7aqezk7AHDRRO5UyhJG/Yf0cKZ/GcGAAAQD6N19GIlPT1dGzZs0AsvvKBnnnlGV199tb71rW/pxhtvHPY9xcXFOuOMMyRJ1113nX7wgx/o9ttv7zv/4osv6n3ve5/S0tIkSe9///v1wgsvyBijK6+8sq8bmZNzdFHG2tpaXXbZZXr44Ye1dOnSGHyniWs8ncKTjDH5/b7+F2PMn4wx3zfG5Iz0XqdyBZLl9XfIFaFT6CuaJf/hStlgcMDxuQty5Q359GjnCVKwS9rx15Efsufp8AIzGX2/Gm093KT01K4BK49KUigB5xRKUklumnLSfMwrBAAAcCi3261zzz1XX/va1/SjH/1IDz/88IjXD/6bdvDXdpgFG60d/u/hrKwsFRcX66WXXhpH5VPDeOYU/q+kCyTJGHO2pG9Juk1SmaS1kq6MdnHHO3cwWZ5gjdyZBUPOeYuKJL9fgZoaeQuOnl+6ZJ6qH9+vHfUeKatY2vJHaeUHIz+gu006sF46+Za+Q+3dAe2ra9O8wg5l91t5VJKUgKuPSuH/Ea+ePe24CoXWWtW1dmt3TavK69vU0NatI23dOtLuV2cgqGDQKhAKyR+0CvX8PyWXMTKm51Xh79tldPSYCR8zPde6+n0tSTJSTqpP/+/dJ8TpuwYAAIi+HTt2yOVyaeHChZKkzZs3a07Pav3DOXDggNavX6/TTjtNDz74oM4888wB588++2zdeOON+uIXvyhrrR555BH96le/ks/n0/ve9z599rOf1fTp09XQ0NDXLfT5fHr00Ud18cUXKz09Xddee21svuEENJ5Q6LbWNvR8frWktdbahyU9bIzZHPXKpgB3MDk8pzDi8NHwikr+gwcHhMJFJSX6m2un0tpDai97t1I3/Vxqq5PScofcQ/uel4Ld0vzz+g5tr2qRtVLItGha0rwBlyfq8FEpPIR03dvVqm/t0vRBW3gkiqYOv/62pUov7q7Ty3vqVdfaNeB8itetaaleJfvc8riMPC6XvG4jl8vI2vBqTNba8O9n0KuVVcgePW8VPtd7Xgr//mZmJubPBgAA4Fi1trbqtttuU2NjozwejxYsWKC1a9eO+J7S0lLdf//9+td//VctXLhQH//4xwecX716tW688UadfHJ43Y2PfvSjWrVqlSTpy1/+ss455xy53W6tWrVK9913X9/70tLS9Pjjj+vCCy9UWlqaLrvssuh+swlqXKHQGOOx1gYknS/pln7norWK6dRhrdyhZHkCww8flaTug4eU2m8Sq9vtUntWk2Z2p+uV7HfqvNBPpE2/ks6MsBXkW7+XUnKkuWf1Hdp6OLwVRWeoZcDKo1LPlhQJ2CmUpJNLwv9C8+q+Br1z+dDOajztqGrRPc/u1hNbqtQdCGlGRpLOXDBdK4uztWBGukpy05SbnqRkrzvepQIAABx31qxZM+zWE5FWHp07d662bRu6Nsfg6z/3uc/pc5/73JBrbrjhBt1www1D7rllyxZJUnZ2tv75z3+Ouf6pYDxh7kFJzxlj6iR1SHpBkowxCyQ1xaC241qgs01u65M7GLlT6CkslIwZsgKpJKUVGqW+nacnKlN13tyzpH/eK532ScntPXpRV4u0/a9S2bWSx9d3+I2KRk1Lc6vV36Ls5OwB903kTuGKoiyl+dx6eU9dwoTCmpZOffPxt/XYm4eV5vPompOKdcWaIi2flZWQczMBAACAYzGe1UfvNMY8JalA0t/t0dmbLoXnFqKf5voaSZIn0ClXxtAVNV0+nzwzZ0YMhYUl01Sz1aMte3ZK77tV+u010uZfS2tuPHrRxl9JgQ5p5TUD3rvpwBEtn+3RJkk5SQMXmrEJOqdQkrxul04uydHLe+rjXYok6dFNh/Sff96qDn9QHz9nvm45e56yU32jvxEAAAA4zox59VFJsta+Yq19xFrb1u/YTmvtxuiXdnyrrq2UJHkCHXJnZUW8xls0S92HhobCpYvDcwG9bUdUkXeOVHSS9Mx/SW09gam7TXrxe1LJOVLx0aGnje3d2lPbpgX54eAXuVOYmKFQkk6fn6u9tW2qauqMWw2BYEhf/fNWfeZ3m7VgRrqe+PRZ+vdLlhAIAQAAMGWNKxRi7OrqwwHOE+yUOzPyPoO+WUXyHzw05PjiefPkd3VpZtCl9XsbpHfeLXU0SH/4sNR4QPrTJ6W2Wum8/zfgfZsqGiVJxbkhSRqyJYUN2YQdPipJp82fLklav7cuLs9v7wroI/e9pgde2qObTp+t33xkjeZOS1LA7w9/dHcr0N0tf3fX0Y+uzvBHZ/iju7Mj/NHR3vcBAAAAJDIWiImRxsZGSely2S6Z5OSI13iLihT4858V6u6Wy3e0E+V2u9Q1rUn5nel6aU+dPnDSKuk935f+9Anpf5aHL7rga1LxyQPut2n/EbmMlJPplyRlJ2UPOJ/Iw0cl6YSCTGWlePXy7nq9b1VRxGsC3d1qaahTS12tOlqa1dnaoo6WFnW2tqiztVX+zg4F/N3yd3Up0N3VF+QC3d0KhYKyoZBCoZBsMBh+DYVkQ0GFgiFZG9JKSSslqVz60a8n/j35UlJ1230PTfxGAAAAQIwQCmOkpWflI+MJDTtk01tUJFmrwOHD8s2dO+Bccr5R0tt5emxHpYKhMrnLrpXyV0j7npNmnSjNPmXI/TYeaNTi/Ex1BMPdx4irjyZuJpTLZXTavOl6eU+9WurrVHegXHUV+1V/8IDqDx5Qc12t2psaI77Xm5SspPR0JaWkyuPzyePzyZeSqtSsaX1fu9xuuVwumZ4Pl8st43LJGqO/vFWtA0c6dNGyQi0tzOpbkWfI72644/2P9V4jyeXxDrkOAAAASCSEwhhpawvPi3P7hk9hvuJwN6y74uCQUFhQkq3arV55uw9o44EjOmlujpS/LPwRQVcgqNf3N+iDJ83Wkc7wcrpZSQPnMoZCSthOYXNdrfZtel3Ldr+mmbu3a+2tLX3n0qblaHrRbM0/cZ4yp+cpIzdPGdNzlZqVreT0DCWnpcvjO7Y5f9ZaffHht/R7Veg7N6/QB04sjta3BAAAABwXCIUx0tnll0uSN234sOItCodCf4TFZk5YPE/PPV6hWaZD696uDofCEWzc36hOf0hnLsjVay1HlOnLlNc1sEtlQ1auBAqFDYcPaftLz2rP66+ppnyPJCk5M1vVvjzlnXqBLr/gFE0vnqOU9KGrt0bL2uf36nevV+iT71hAIAQAAIAjsdBMjPi7wou9JGemDXuNZ8YMGa834rYUS3oWm5nr8uipt2tGfd5Lu+vkdhmdMi9HRzqPDBk6KiXG5vUBv19bn3tKv/3PL+iXn/1XvfLw7+RNTtLZH/qwbvzuT3Tr2l+pfOWVesF3gopKl8U0EL68u07fenK73rWiQJ+7cFHMngMAAIDhlZeXa9myyKPh3G63ysrKtGzZMl111VVqbw8v4nfnnXdq6dKlWrFihcrKyvTqq69OZslTDp3CGPF3S0k2qNTs7GGvMS6XvIWF6o6wAqnH7VbntEbldKZqd02rdlS1aHH+8AHp2Z01KivOVkayV0e6jmha0tBQGApZeeMUCrs72vXmuif1+l8eVduRBk0rKNRZ196oE84+T+nTBnZBzy+doZ8+t1dN7X5lpcZmTl5je7c+99AbKslN091XrkyoDioAAADCUlJStHnzZknShz70If30pz/Vaaedpscff1wbN25UUlKS6urq1N3dHd9Cj3OEwhgJBV1yBzqUXJA94nXeoiL5Dw0NhZKUnC+lvD1d3iy//rjpoL50aWnE6yoa2rXlULO+dOkSSdKRziMqTC8ccl08OoXBQEBv/OMJrf/Db9TZ2qLZy1bokls/qznLy4ZdgOf80pn68TN79OzOGl1WNivqNVlr9R+PvKX6ti79/IYzlOJzR/0ZAAAAx5tn7lurmv17o3rPGXPm6R033jLqdcFgUDfffLNefvllzZo1S3/605+UkpIy4JqzzjpLb775pubOnavc3FwlJSVJknJzc6NasxMxfDRGQn6PvIHh9yjs5SnIl7+qMuK5grnT5An5dHZxQI9uOqRgyEa87sktVZKkS5cVSJIaOxuH7FEohTuFkxkKy9/cpPs//0k9c9//asbcebr2m/+tq75yl+auWDVsIJSksqJs5ab79I9t1TGp6w8bDuqvb1Xpcxcu1rJZWaO/AQAAADG1a9cufeITn9DWrVuVnZ2thx9+eMD5QCCgJ554QsuXL9dFF12kiooKLVq0SLfeequee+65OFU9ddApjBFXwBveuD5r5FDoLShQsLZuyF6FknTC4hI9/5dDWpDs11PNXXphV63OXTxjwDXWWj26+ZCWzcrU7OmpstaqoathyB6F0uQtNNPZ1qrnfvULbXnmH5pWUKjL//0OzVt90ohBsD+Xy+iipfl6ZOMhtXUFlJYUvf9M99e36at/3qpTSnJ0y9nzonZfAACA491YOnqxUlJSorKyMknSmjVrVF5eLknq6OjoO37WWWfppptuks/n04YNG/TCCy/omWee0dVXX61vfetbuvHGG+NS+1RAKIwRVyBJ7kCnXJlD5/b15y0ID/MMVFXJN3v2gHOl8xfoKfdeeZq6NDMzR2uf3zskFG480Kith5t15/vCk3Pb/G0KhAKRO4XB2HcKD769RX/5wf+ntiNHdPJlV+q0K689pu0iLi+bpd+8ekD/2Faty1dFZwhpIBjSZ363WS6X0XevLpObeYQAAAAJoXcoqBReXKajo0PSwDmF/bndbp177rk699xztXz5ct1///2Ewglg+GiM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
"
>
</div>

</div>

</div>

</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[815]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">sld_si</span> <span class="o">=</span> <span class="n">s_d2oblock</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">()</span>
<span class="n">sld_simucinsd2o</span> <span class="o">=</span> <span class="n">s_d2omucins</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">()</span>
<span class="n">sld_simucinsh2o</span> <span class="o">=</span> <span class="n">s_h2omucins</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">()</span>
<span class="n">sld_siconfimucinsd2oP1</span> <span class="o">=</span> <span class="n">s_d2omucinsconfined1</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">()</span>
<span class="n">sld_siconfimucinsd2oP2</span> <span class="o">=</span> <span class="n">s_d2omucinsconfined2</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">()</span>
<span class="n">sld_hps</span> <span class="o">=</span> <span class="n">s_hPS</span><span class="o">.</span><span class="n">sld_profile</span><span class="p">()</span>

<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sld_si.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sld_si</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sld_simucinsd2o.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sld_simucinsd2o</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sld_simucinsh2o.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sld_simucinsh2o</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sld_siconfimucinsd2oP1.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sld_si</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sld_siconfimucinsd2oP2.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sld_simucinsd2o</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
<span class="n">np</span><span class="o">.</span><span class="n">savetxt</span><span class="p">(</span><span class="s1">&#39;sld_hps.txt&#39;</span><span class="p">,</span><span class="n">np</span><span class="o">.</span><span class="n">transpose</span><span class="p">(</span><span class="n">sld_simucinsh2o</span><span class="p">),</span><span class="n">fmt</span><span class="o">=</span><span class="s1">&#39;</span><span class="si">%.18e</span><span class="s1">&#39;</span><span class="p">,</span> <span class="n">delimiter</span><span class="o">=</span><span class="s1">&#39; &#39;</span><span class="p">)</span>
</pre></div>

     </div>
</div>
</div>
</div>

</div><div class="jp-Cell jp-CodeCell jp-Notebook-cell jp-mod-noOutputs  ">
<div class="jp-Cell-inputWrapper">
<div class="jp-InputArea jp-Cell-inputArea">
<div class="jp-InputPrompt jp-InputArea-prompt">In&nbsp;[&nbsp;]:</div>
<div class="jp-CodeMirrorEditor jp-Editor jp-InputArea-editor" data-type="inline">
     <div class="CodeMirror cm-s-jupyter">
<div class=" highlight hl-ipython3"><pre><span></span> 
</pre></div>

     </div>
</div>
</div>
</div>

</div>
</body>







</html>